#!/usr/bin/env python3
"""Fetch (or generate) NACA 4-digit airfoil coordinate files.

Primary source: UIUC Airfoil Database (public domain). If the host is
unreachable, coordinates are generated analytically from the NACA 4-digit
thickness/camber equations with cosine (half-integer) spacing, which is
equivalent to the database content for these profiles.

Output format (Selig/LEDNICER single-loop, trailing edge first):
    x  y        (upper surface TE -> LE, then lower surface LE -> TE)
"""
import os
import sys

import numpy as np

AIRFOILS = {
    "naca0012": "https://m-selig.ae.illinois.edu/ads/coord/naca0012.dat",
    "naca2412": "https://m-selig.ae.illinois.edu/ads/coord/naca2412.dat",
    "naca4412": "https://m-selig.ae.illinois.edu/ads/coord/naca4412.dat",
}

N_POINTS = 320  # per surface; dense enough for any downstream re-sampling


def naca4_digit_string(code: str):
    m, p, tt = int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0
    return m, p, tt


def thickness_distribution(x, t):
    return 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                    + 0.2843 * x**3 - 0.1015 * x**4)  # open TE (NACA TR-824)


def camber_line(x, m, p):
    if m == 0 or p == 0:
        return np.zeros_like(x), np.zeros_like(x)
    yc = np.where(x < p,
                  m / p**2 * (2 * p * x - x**2),
                  m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * x - x**2))
    dyc = np.where(x < p,
                   2 * m / p**2 * (p - x),
                   2 * m / (1 - p)**2 * (p - x))
    return yc, dyc


def naca4_coordinates(code: str, n=N_POINTS, closed_te=False):
    """Return (x, y) single-loop coordinates, TE -> upper -> LE -> lower -> TE."""
    m, p, t = naca4_digit_string(code)
    beta = np.linspace(0.0, np.pi, n)
    x = 0.5 * (1.0 - np.cos(beta))          # cosine clustering at LE and TE
    yt = thickness_distribution(x, t)
    if closed_te:
        yt[-1] = 0.0
    yc, dyc = camber_line(x, m, p)
    theta = np.arctan(dyc)
    xu, yu = x - yt * np.sin(theta), yc + yt * np.cos(theta)
    xl, yl = x + yt * np.sin(theta), yc - yt * np.cos(theta)
    x_loop = np.concatenate([xu[::-1], xl[1:-1]])
    y_loop = np.concatenate([yu[::-1], yl[1:-1]])
    return x_loop, y_loop


def normalise(x, y):
    """Place leading edge at origin, trailing edge at (1, 0), chord = 1."""
    i_le = np.argmin(x)
    c = x.max() - x.min()
    x = (x - x[i_le]) / c
    y = (y - y[i_le]) / c
    # rotate so that TE sits on y = 0
    x_te, y_te = x[np.argmax(x)], y[np.argmax(x)]
    ang = np.arctan2(y_te, x_te)
    xr = x * np.cos(-ang) - y * np.sin(-ang)
    yr = x * np.sin(-ang) + y * np.cos(-ang)
    return xr - xr.min() * 0, yr  # LE already at ~origin


def write_dat(path, x, y, header):
    with open(path, "w") as f:
        f.write(header + "\n")
        for xi, yi in zip(x, y):
            f.write(f"{xi:12.7f}  {yi:12.7f}\n")


def fetch(name, url, out_dir):
    try:
        import requests
        r = requests.get(url, timeout=20)
        if r.status_code == 200 and len(r.text) > 500:
            open(os.path.join(out_dir, f"{name}.dat"), "w").write(r.text)
            print(f"[fetch] {name}: downloaded from UIUC ({len(r.text)} B)")
            return True
        print(f"[fetch] {name}: HTTP {r.status_code} from UIUC")
    except Exception as exc:
        print(f"[fetch] {name}: {type(exc).__name__}: {exc}")
    return False


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "airfoils"
    os.makedirs(out_dir, exist_ok=True)
    for name, url in AIRFOILS.items():
        code = name.replace("naca", "")
        if fetch(name, url, out_dir):
            continue
        x, y = naca4_coordinates(code)
        write_dat(os.path.join(out_dir, f"{name}.dat"), x, y,
                  f"NACA {code.upper()} generated analytically (NACA TR-824 equations)")
        print(f"[gen]   {name}: generated analytically ({len(x)} pts)")


if __name__ == "__main__":
    main()
