#!/bin/bash
# URANS (pimpleFoam) runner for the low-Reynolds-number LSB campaign.
# Two-stage integration: stage 1 settles the gamma-Re_thetat model with a
# fixed tiny dt (cold/warm start survival), stage 2 integrates with
# Co-limited adaptive dt to capture LSB/shedding dynamics. Runs 4-rank
# parallel on the 50c medium mesh; initial condition copied from the
# converged steady sibling (identical grid connectivity).
# Usage: bash scripts/09_run_urans.sh <manifest>   (one case dir per line)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export MAMBA_ROOT_PREFIX=/workspace/tools/mamba/root
set +u
source /workspace/tools/env/etc/bashrc >/dev/null 2>&1
source /workspace/tools/env/etc/bashrc >/dev/null 2>&1
set -u

MANIFEST="${1:-$ROOT/run_manifest_urans.txt}"

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

run_urans_case() {
    local case_dir="$1"
    local base="$(basename "$case_dir")"
    local log="$ROOT/logs/${base}.log"
    {
        echo "=== $(date) : $case_dir ==="
        if [ ! -f "$case_dir/constant/polyMesh/points" ]; then
            gmshToFoam "$(cat "$case_dir/mesh_source.txt")" \
                -case "$case_dir" < /dev/null > /dev/null 2>&1 || { echo "gmshToFoam FAILED"; exit 1; }
            fix_empty_patches "$case_dir"
        fi
        checkMesh -case "$case_dir" < /dev/null 2>&1 | grep -E "Mesh OK|non-orthogonality|skewness|Failed" || true
        local sib="$ROOT/openfoam_cases/$(echo "$base" | sed 's/_urans$//')"
        if [ ! -d "$sib" ]; then
            sib="$ROOT/openfoam_cases/$(echo "$base" | sed 's/_urans$//; s/aoa[0-9]*/aoa10/')"
            echo "no exact steady sibling; falling back to $sib"
        fi
        local last="$(ls -d "$sib"/[1-9]* 2>/dev/null | sort -n | tail -1)"
        if [ -n "$last" ] && [ -f "$last/U" ]; then
            for f in U p k omega nut ReThetat gammaInt; do
                cp "$last/$f" "$case_dir/0/$f" 2>/dev/null || true
            done
            echo "warm-started by direct field copy from $sib/$last"
        fi
        cp "$case_dir/system/controlDict.stage1" "$case_dir/system/controlDict"
        if ls -d "$case_dir"/0.005 >/dev/null 2>&1; then
            echo "stage1 already complete (0.005 exists); skipping"
        else
            pimpleFoam -case "$case_dir" < /dev/null > /dev/null 2>&1
            echo "stage1 rc=$? $(date)"
        fi

        if ls -d "$case_dir"/0.005 >/dev/null 2>&1; then
            cp "$case_dir/system/controlDict.stage2" "$case_dir/system/controlDict"
            pimpleFoam -case "$case_dir" < /dev/null > /dev/null 2>&1
            echo "stage2 rc=$? $(date)"
        else
            echo "stage2 SKIPPED: no 0.005 field directory from stage1"
        fi
        echo "=== $(date) : finished ==="
    } > "$log" 2>&1
    grep -E "^===|Mesh OK|Failed|warm-started|stage" "$log" | tail -6
}

export -f run_urans_case fix_empty_patches
export ROOT

mapfile -t CASES < "$MANIFEST"
for c in "${CASES[@]}"; do run_urans_case "$c"; done
echo "All URANS cases processed."
