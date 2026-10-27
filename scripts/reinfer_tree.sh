#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/generated/reinferred_tree"
mkdir -p "$OUT"
# Requires MAFFT, trimAl, and IQ-TREE 3; see README for version limitations.
# A separate output directory protects the archived manuscript tree.
mafft --localpair --maxiterate 1000 "$ROOT/data/figure_E/TerL_all.faa" > "$OUT/TerL_aligned.faa"
trimal -in "$OUT/TerL_aligned.faa" -out "$OUT/TerL_trimmed.faa" -automated1
iqtree3 -s "$OUT/TerL_trimmed.faa" -m MFP -B 1000 --alrt 1000 -seed 886838 -T 1 --prefix "$OUT/TerL_reinferred"
