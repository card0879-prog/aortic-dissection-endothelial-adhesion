# Code — Endothelial adhesion molecules in human aortic dissection (multicohort single-cell)

Reproduces the analyses in the manuscript (patient-level pseudobulk + random-effects meta-analysis
of SELE/SELP/VCAM1/ICAM1 in aortic dissection endothelium). Computational only; no wet-lab data.

## Data access (not redistributed here)
Public scRNA-seq: GSE155468, GSE189795, GSE213740, GSE216860, GSE222318 (NCBI GEO).
In-house scRNA-seq: CNSA/CNGBdb accession CNP0007139.
Analyses start from an endothelial subset object `ec_subset.rds` (VEC+LEC, with meta columns
celltype, subtype_2, Source, Group, orig.ident), derived from the integrated atlas (Tao et al.,
Adv Sci 2026, doi:10.1002/advs.75509).

## Run order (R; cluster via `Rscript`)
`00_sessionInfo.R` records the R session information. The environment used for the analyses has already been documented in `env/sessionInfo.txt` and `env/r_packages.txt`.
01_pseudobulk_meta.R   AUTHORITATIVE (reproduces reported tables): cohort counts -> S1; per-cohort voom -> S2A; metafor REML+KNHA meta -> S2B; leave-one-out (HK+conventional) -> S4; within-subtype 12 tests -> S5; design rank check + per-source counts -> S5b
02_edgeR_DESeq2.R      multi-method sensitivity edgeR-QLF + DESeq2 per cohort -> S3
03_KLF2_sensitivity.R  per-cohort adhesion AD coef with/without KLF2 covariate (wide) -> TableS6_KLF2_sensitivity.csv
04_detection.R            per-cohort per-gene detection fraction + mean logCPM -> TableS7
05_fig1_coords.R          export UMAP coords for Fig 1 -> fig1_coords.csv

## Figures (Python)
`python/meta.py` provides an approximate independent cross-check of pooled effects; it is not the authoritative source for the reported meta-analysis statistics. `python/make_figures.py` generates Figures 1–4 and Supplementary Figure S1 from the supplied summary data. All reported pooled estimates, confidence intervals, heterogeneity statistics, and p-values are taken from the R `metafor` analysis.

## Environment
R 4.2.2 with Seurat 5.0.2, limma 3.54.2, edgeR 3.40.2, DESeq2 1.38.3, metafor 4.0.0, and Matrix 1.6-4. Full session details are recorded in `env/sessionInfo.txt`. 
Python dependencies are listed in `env/requirements.txt`.

## Notes
- Patient-level aggregation collapses in-house proximal/medial/distal segments of one patient
  (orig.ident "ID-A/B/C") into a single patient to avoid pseudoreplication.
- Per-cohort voom order (01): filterByExpr (force-retaining the four genes) THEN calcNormFactors; SE = stdev.unscaled[,2]*sqrt(s2.post). See TableS7 for detection rates.
- Summary outputs (TableS1–S7) are the non-sensitive results underlying all figures.

## License
Code: MIT (see LICENSE). Summary tables (TableS1–S7) and figures: CC-BY 4.0.

## Citation
See CITATION.cff. Please cite the article and this repository.

## Suggested GitHub "About"
Description: Patient-level multicohort single-cell meta-analysis of endothelial
adhesion-molecule expression in human aortic dissection (R + Python; reproduces
all figures and supplementary tables).
Topics: single-cell, aortic-dissection, meta-analysis, bioinformatics, rna-seq,
endothelium, reproducible-research

## Reproducibility & data tiers
Three tiers of inputs:
1. **Public raw scRNA-seq** (GEO: GSE155468 [Li 2020], GSE189795 [He 2022], GSE213740 [Zhang 2023], GSE222318 [Zhao 2025], GSE216860 [accession record]) — download from GEO and preprocess following the integrated-atlas pipeline (Tao et al. 2026, doi:10.1002/advs.75509: Seurat QC, Harmony integration, annotation) to obtain the endothelial subset.
2. **Restricted intermediate `ec_subset.rds`** (patient-level endothelial subset with the annotations used here) is NOT redistributed: it derives in part from in-house scRNA-seq governed by the original project's data-sharing terms (available via CNSA/CNGBdb CNP0007139 on reasonable request). R scripts 01–05 start from this object.
3. **Non-sensitive summary outputs** — `summary_data/` (TableS1–S7 + fig1_coords.csv) — ARE included here.

**What reproduces from what**
- From the restricted `ec_subset.rds`: `R/01–05` regenerate every TableS* from counts (pass the data dir as arg1 or set `$AD_ATLAS_DIR`; no hard-coded paths).
- From the public `summary_data/` ALONE (no restricted data): all figures and the independent statistical cross-check —
  - `python python/make_figures.py --dir summary_data --out figures`  (Fig1–4, FigS1)
  - `python python/meta.py summary_data/TableS2A_percohort_voom.csv`  (independent pooled-effect check)

**Authoritative statistics**: all pooled effects, CIs, I² and p-values come from **R metafor** (`R/01` → TableS2B/S4). `python/meta.py` is an independent approximate cross-check only; its Q-based I² differs from metafor's REML I² and must not be used.

**Statistical precision**: `R/01` passes full-precision per-cohort log2FC/SE to metafor (rounding only on display), matching the reported tables; a perturbation check confirms 3-dp rounding does not change any reported significance (ICAM1 p varies <0.0002).

**Within-subtype design check**: `R/01` performs a design-matrix rank check for `~Source + Group` (falling back to `~Group` and flagging if rank-deficient) and writes per-subtype, per-source AD/Ctrl patient counts to `TableS5b_subtype_source_counts.csv`.

**Software environment:** The R package versions used for the analyses are recorded in `env/r_packages.txt`, with full session information in `env/sessionInfo.txt`.
