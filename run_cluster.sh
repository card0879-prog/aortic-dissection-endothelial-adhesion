#!/usr/bin/env bash
# Regenerate all supplementary tables from ec_subset.rds.
# Usage:  AD_ATLAS_DIR=/path/with/ec_subset.rds  bash run_cluster.sh
set -e
DIR="${AD_ATLAS_DIR:-.}"
RS="${RSCRIPT:-Rscript}"
"$RS" R/00_sessionInfo.R
for s in 01_pseudobulk_meta 02_edgeR_DESeq2 03_KLF2_sensitivity 04_detection 05_fig1_coords; do
  "$RS" "R/${s}.R" "$DIR"
done
echo "ALL DONE — tables written to working dir; paste sessionInfo.txt into env/r_packages.txt"
