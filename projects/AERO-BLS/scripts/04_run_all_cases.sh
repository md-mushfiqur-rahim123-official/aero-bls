#!/bin/bash
# Batch runner: converts mesh, checks quality, runs simpleFoam for every case
# registered in run_manifest.txt (one case directory per line).
# Usage: bash scripts/04_run_all_cases.sh [n_parallel]
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export MAMBA_ROOT_PREFIX=/workspace/tools/mamba/root
set +u
source /workspace/tools/env/etc/bashrc >/dev/null 2>&1
set -u

NPAR="${1:-1}"
MANIFEST="${MANIFEST:-$ROOT/run_manifest.txt}"

fix_empty_patches() {
    local bf="$1/constant/polyMesh/boundary"
    python3 - "$bf" <<'PYEOF'
import sys, re
path = sys.argv[1]
text = open(path).read()
for name in ("front", "back"):
    text = re.sub(r"(%s\s*\{[^}]*?type\s+)patch;" % name, r"\1empty;", text)
text = re.sub(r"(airfoil\s*\{[^}]*?type\s+)patch;", r"\1wall;", text)
open(path, "w").write(text)
PYEOF
}

run_case() {
    local case_dir="$1"
    local log="$ROOT/logs/$(basename "$case_dir").log"
    {
        echo "=== $(date) : $case_dir ==="
        if [ ! -d "$case_dir/constant/polyMesh" ] || [ ! -f "$case_dir/constant/polyMesh/points" ]; then
            gmshToFoam "$(cat "$case_dir/mesh_source.txt")" \
                -case "$case_dir" > /dev/null 2>&1 || { echo "gmshToFoam FAILED"; exit 1; }
            fix_empty_patches "$case_dir"
        fi
        checkMesh -case "$case_dir" 2>&1 | grep -E "Mesh OK|non-orthogonality|skewness|Failed" || true
        simpleFoam -case "$case_dir" > /dev/null 2>&1
        local rc=$?
        echo "=== $(date) : finished rc=$rc ==="
    } > "$log" 2>&1
    grep -E "^===|Mesh OK|Failed" "$log" | tail -3
}

export -f run_case fix_empty_patches
export ROOT

mapfile -t CASES < "$MANIFEST"
if [ "$NPAR" -le 1 ]; then
    for c in "${CASES[@]}"; do run_case "$c"; done
else
    printf "%s\n" "${CASES[@]}" | xargs -P "$NPAR" -I{} bash -c 'run_case "$@"' _ {}
fi
echo "All cases processed."
