#!/usr/bin/env bash
# One-time environment setup for the AERO-BLS pipeline.
# Targets WSL2 Ubuntu (as documented in the project brief) and plain Debian/Ubuntu.
# Installs: Python scientific stack, Gmsh (pip), OpenFOAM + tectonic (micromamba/conda-forge).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS="${TOOLS_DIR:-$ROOT/../tools}"
ENV="$TOOLS/env"
MAMBA_ROOT_PREFIX="$TOOLS/mamba/root"
export MAMBA_ROOT_PREFIX

echo "== AERO-BLS setup =="
echo "Project root: $ROOT"
echo "Toolchain   : $ENV"

# ---- 1. System packages (needs sudo; on managed systems use the micromamba path below) ----
if command -v sudo >/dev/null 2>&1 && command -v apt-get >/dev/null 2>&1; then
  sudo apt-get update -qq || true
  sudo apt-get install -y -qq python3 python3-pip python3-venv curl bzip2 || true
fi

# ---- 2. Python stack ----
python3 -m pip install --user --quiet numpy scipy matplotlib pyvista meshio requests || true

# ---- 3. micromamba (user-space conda) ----
if [ ! -x "$TOOLS/mamba/bin/micromamba" ]; then
  mkdir -p "$TOOLS/mamba/bin"
  echo "-- downloading micromamba --"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest -o /tmp/mm.tar.bz2
  tar -xjf /tmp/mm.tar.bz2 -C "$TOOLS/mamba/bin" bin/micromamba --strip-components=1
fi
MM="$TOOLS/mamba/bin/micromamba"

# ---- 4. OpenFOAM (ESI v2412) + tectonic LaTeX ----
if [ ! -x "$ENV/bin/simpleFoam" ]; then
  echo "-- installing openfoam + tectonic into $ENV (this takes a few minutes) --"
  "$MM" create -y -p "$ENV" -c conda-forge openfoam=2412 tectonic
fi

# ---- 5. Gmsh python module ----
python3 -m pip install --user --quiet gmsh || true

cat <<EOF

== Setup complete ==
Activate the toolchain in every new shell:

  export MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX
  set +u; source $ENV/etc/bashrc; set -u   # if your shell has 'set -u' active

Verify:
  $ENV/bin/simpleFoam -help | head -1
  $ENV/bin/tectonic --version
  $ENV/bin/blockMesh -help | head -1
EOF
