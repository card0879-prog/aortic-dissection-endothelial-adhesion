# Capture exact package versions for reproducibility.
suppressMessages({library(Seurat);library(limma);library(edgeR);library(DESeq2);library(metafor);library(Matrix)})
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("wrote sessionInfo.txt\n")
