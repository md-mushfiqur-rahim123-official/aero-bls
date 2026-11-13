#!/usr/bin/env python3
"""Generate a structured O-grid (C-type topology) around a 2D airfoil.

The grid is fully structured: N_theta rays emanating from an interior
centre point, each carrying n_layers cells stretched geometrically from
the wall (first-cell height ~1e-5 c => y+ < 1 at Re = 2e6) to a circular
farfield of R chord lengths. The mesh is written in Gmsh ASCII format
(v2.2) as a single-cell-thick hexahedral grid, ready for `gmshToFoam`.

Optional 2D surrogates of 3D passive devices (documented limitations):
  --serration  sawtooth trailing-edge notches (2D surrogate of TE serrations)
  --tubercle   sinusoidal leading-edge thickness modulation (2D surrogate)
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load_airfoil(path):
    xs, ys = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                try:
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))
                except ValueError:
                    continue
    x, y = np.array(xs), np.array(ys)
    # normalise: LE at origin, chord = 1, chord line on y = 0
    i_le = np.argmin(x)
    x, y = x - x[i_le], y - y[i_le]
    c = x.max()
    x, y = x / c, y / c
    ang = np.arctan2(y[np.argmax(x)], x[np.argmax(x)])
    ca, sa = np.cos(-ang), np.sin(-ang)
    return x * ca - y * sa, x * sa + y * ca


def naca4_surface_functions(code):
    m, p, t = int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0

    def yt(x):
        return 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                        + 0.2843 * x**3 - 0.1015 * x**4)

    def yc(x):
        if m == 0 or p == 0:
            return np.zeros_like(np.asarray(x, dtype=float))
        return np.where(x < p,
                        m / p**2 * (2 * p * x - x**2),
                        m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * x - x**2))

    def dyc(x):
        if m == 0 or p == 0:
            return np.zeros_like(np.asarray(x, dtype=float))
        return np.where(x < p,
                        2 * m / p**2 * (p - x),
                        2 * m / (1 - p)**2 * (p - x))
    return yt, yc, dyc


def apply_serration(x, y, x0=0.85, n_teeth=5, depth_frac=0.6):
    """Sawtooth notches cut into the trailing-edge region (2D surrogate)."""
    yt, yc, dyc = None, None, None
    half_th = np.zeros_like(x)
    # local half-thickness estimated from the loop itself
    upper = y >= 0
    for s in (upper, ~upper):
        xs = x[s]
        order = np.argsort(xs)
        half_th[s] = np.interp(x[s], xs[order], np.abs(y[s][order]))
    lam = (1.0 - x0) / n_teeth
    frac = np.mod((x - x0) / lam, 1.0)
    tri = np.where(x >= x0, 1.0 - np.abs(2.0 * frac - 1.0), 0.0)
    tri *= (x >= x0) & (x <= 1.0)
    depth = depth_frac * half_th
    y = np.where(upper, y - depth * tri, y + depth * tri)
    return x, y


def apply_tubercle(x, y, amplitude=0.004, n_bumps=3, x_end=0.30):
    """Sinusoidal leading-edge thickness modulation (2D surrogate)."""
    w = np.clip((x_end - x) / (0.05 * x_end), 0.0, 1.0)
    w = w**2 * (3 - 2 * w)
    pert = amplitude * np.sin(2 * np.pi * n_bumps * x / x_end) * w
    return x, y + pert


def resample_by_angle(x, y, n_theta, centre=(0.4, 0.0)):
    """Resample the closed loop uniformly (cosine-clustered) in polar angle
    about an interior centre; requires a star-shaped profile."""
    cx, cy = centre
    th = np.mod(np.arctan2(y - cy, x - cx), 2.0 * np.pi)
    # enforce a strictly increasing branch along the loop (single 2pi wrap)
    d = np.diff(th)
    d[d < -np.pi] += 2.0 * np.pi
    d[d > np.pi] -= 2.0 * np.pi
    th = th[0] + np.concatenate([[0.0], np.cumsum(d)])
    # cosine clustering on each half: dense at TE (both ends of the loop)
    # and at the leading edge (segment junction)
    i_le = int(np.argmin(x))
    th_le = th[i_le]
    n1 = n_theta // 2
    n2 = n_theta - n1
    t1 = 0.5 * (1.0 - np.cos(np.pi * np.arange(n1 + 1) / n1))
    th_a = th[0] + (th_le - th[0]) * t1            # TE-upper -> LE
    t2 = 0.5 * (1.0 - np.cos(np.pi * np.arange(1, n2 + 1) / n2))
    th_b = th_le + (th[-1] - th_le) * t2           # LE -> TE-lower
    th_t = np.concatenate([th_a, th_b])
    xs = np.interp(th_t, th, x)
    ys = np.interp(th_t, th, y)
    return xs, ys, th_t


def build_grid(surf_x, surf_y, n_theta, n_layers, r_far, first_cell, growth,
               centre=(0.5, 0.0), j_blend=30):
    cx, cy = centre
    th = np.arctan2(surf_y - cy, surf_x - cx)
    far_x = cx + r_far * np.cos(th)
    far_y = cy + r_far * np.sin(th)
    seg_len = np.hypot(far_x - surf_x, far_y - surf_y)
    # geometric stretching: sum d_k = seg_len, d_k = d0 * growth**k
    d0 = seg_len * (growth - 1.0) / (growth**n_layers - 1.0)
    t = (growth ** np.arange(n_layers + 1) - 1.0) / (growth**n_layers - 1.0)
    Xp = surf_x[None, :] + (far_x - surf_x)[None, :] * t[:, None]
    Yp = surf_y[None, :] + (far_y - surf_y)[None, :] * t[:, None]

    # near-wall layers follow surface normals (orthogonal BL), blended into
    # the polar grid over j_blend layers (j_blend=0 disables blending)
    tx = np.gradient(surf_x)
    ty = np.gradient(surf_y)
    ln = np.hypot(tx, ty)
    nx = ty / ln
    ny = -tx / ln
    flip = (nx * (surf_x - cx) + ny * (surf_y - cy)) < 0
    nx = np.where(flip, -nx, nx)
    ny = np.where(flip, -ny, ny)
    s = t[:, None] * seg_len[None, :]
    Xn = surf_x[None, :] + nx[None, :] * s
    Yn = surf_y[None, :] + ny[None, :] * s
    if j_blend > 0:
        w = np.clip(np.arange(n_layers + 1) / float(j_blend), 0.0, 1.0) ** 2
        w[0] = 0.0
        X = (1 - w)[:, None] * Xn + w[:, None] * Xp
        Y = (1 - w)[:, None] * Yn + w[:, None] * Yp
    else:
        X, Y = Xp, Yp
    return X, Y, d0


def write_gmsh_v22(path, X, Y, thickness, phys):
    """Write a single-layer hexahedral grid as Gmsh ASCII v2.2.

    phys: dict with physical tags {airfoil:1, farfield:2, back:3, front:4, vol:5}
    """
    n_layers1, n_theta = X.shape           # j = 0..n_layers (inclusive)
    nz = 2
    node_ids = np.arange(1, n_theta * nz * n_layers1 + 1).reshape(nz, n_layers1, n_theta)

    def nid(i, j, k):
        return node_ids[k, j, i % n_theta]

    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n5\n")
        f.write('1 1 "airfoil"\n1 2 "farfield"\n2 3 "back"\n2 4 "front"\n3 5 "internalMesh"\n')
        f.write("$EndPhysicalNames\n")
        f.write(f"$Nodes\n{node_ids.size}\n")
        for k, z in enumerate((0.0, thickness)):
            for j in range(n_layers1):
                for i in range(n_theta):
                    f.write(f"{nid(i, j, k)} {X[j, i]:.9f} {Y[j, i]:.9f} {z:.6f}\n")
        f.write("$EndNodes\n")

        n_wall = n_theta
        n_far = n_theta
        n_cap = n_theta * (n_layers1 - 1)
        n_hex = n_theta * (n_layers1 - 1)
        n_elem = n_wall + n_far + 2 * n_cap + n_hex
        f.write(f"$Elements\n{n_elem}\n")
        e = 0
        # wall & farfield side quads (type 3), outward-facing
        for i in range(n_theta):
            e += 1
            f.write(f"{e} 3 2 1 1 {nid(i,0,0)} {nid(i,0,1)} {nid(i+1,0,1)} {nid(i+1,0,0)}\n")
        for i in range(n_theta):
            e += 1
            f.write(f"{e} 3 2 2 2 {nid(i,-1,0)} {nid(i+1,-1,0)} {nid(i+1,-1,1)} {nid(i,-1,1)}\n")
        # back (z=0) and front (z=T) caps
        for i in range(n_theta):
            for j in range(n_layers1 - 1):
                quad = (nid(i, j, 0), nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j, 0))
                e += 1
                f.write(f"{e} 3 2 3 3 " + " ".join(map(str, quad)) + "\n")
        for i in range(n_theta):
            for j in range(n_layers1 - 1):
                quad = (nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1), nid(i, j + 1, 1))
                e += 1
                f.write(f"{e} 3 2 4 4 " + " ".join(map(str, quad)) + "\n")
        # hexahedra (type 5)
        for i in range(n_theta):
            for j in range(n_layers1 - 1):
                hexn = (nid(i, j, 0), nid(i + 1, j, 0), nid(i + 1, j + 1, 0), nid(i, j + 1, 0),
                        nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1), nid(i, j + 1, 1))
                e += 1
                f.write(f"{e} 5 2 5 5 " + " ".join(map(str, hexn)) + "\n")
        f.write("$EndElements\n")
    return n_hex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--airfoil", required=True, help="path to .dat or nacaXXXX code")
    ap.add_argument("--output", required=True, help="output .msh path")
    ap.add_argument("--density", default="medium", choices=["coarse", "medium", "fine"])
    ap.add_argument("--serration", action="store_true")
    ap.add_argument("--tubercle", action="store_true")
    ap.add_argument("--farfield", type=float, default=20.0)
    ap.add_argument("--first-cell", type=float, default=1e-5)
    ap.add_argument("--centre", type=float, default=0.5)
    ap.add_argument("--centre-y", type=float, default=0.0)
    ap.add_argument("--j-blend", type=int, default=0)
    args = ap.parse_args()

    settings = {"coarse": (160, 100, 1.095),
                "medium": (240, 140, 1.09),
                "fine": (320, 160, 1.085)}
    n_theta, n_layers, growth = settings[args.density]

    name = os.path.basename(args.airfoil)
    if os.path.exists(args.airfoil):
        x, y = load_airfoil(args.airfoil)
    else:
        code = os.path.splitext(name)[0].replace("naca", "")
        yt, yc, dyc = naca4_surface_functions(code)
        beta = np.linspace(0.0, np.pi, 400)
        xs = 0.5 * (1.0 - np.cos(beta))
        h = yt(xs)
        th_ = np.arctan(dyc(xs))
        x = np.concatenate([(xs - h * np.sin(th_))[::-1], (xs + h * np.sin(th_))[1:-1]])
        y = np.concatenate([(yc(xs) + h * np.cos(th_))[::-1], (yc(xs) - h * np.cos(th_))[1:-1]])
    if args.serration:
        x, y = apply_serration(x, y)
    if args.tubercle:
        x, y = apply_tubercle(x, y)

    sx, sy, th = resample_by_angle(x, y, n_theta, centre=(args.centre, args.centre_y))
    X, Y, d0 = build_grid(sx, sy, n_theta, n_layers, args.farfield,
                          args.first_cell, growth, centre=(args.centre, args.centre_y),
                          j_blend=args.j_blend)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    n_hex = write_gmsh_v22(args.output, X, Y, 0.1, None)
    meta = {"airfoil": name, "density": args.density, "n_theta": n_theta,
            "n_layers": n_layers, "cells": int(n_hex),
            "first_cell_height": float(d0.min()), "growth": growth,
            "farfield_radius": args.farfield,
            "serration": args.serration, "tubercle": args.tubercle}
    with open(os.path.splitext(args.output)[0] + "_meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[mesh] {args.output}: {n_hex} cells "
          f"({n_theta} x {n_layers}), dy1={d0.min():.2e}")


if __name__ == "__main__":
    main()
