# Export UMAP coordinates + annotations for Figure 1
suppressMessages(library(Seurat))
args<-commandArgs(trailingOnly=TRUE); setwd(if(length(args)>=1) args[[1]] else Sys.getenv("AD_ATLAS_DIR","."))
ec<-readRDS("ec_subset.rds")
rn<-names(ec@reductions); um<-rn[grepl("umap",tolower(rn))][1]; if(is.na(um)) um<-tail(rn,1)
emb<-Embeddings(ec,um); md<-ec@meta.data
df<-data.frame(UMAP1=emb[,1],UMAP2=emb[,2],subtype=as.character(md$subtype_2),
  celltype=as.character(md$celltype),Group=as.character(md$Group),Source=as.character(md$Source))
set.seed(1); df<-df[sample(nrow(df),min(30000,nrow(df))),]
write.csv(df,"fig1_coords.csv",row.names=FALSE); cat("DONE 05\n")
