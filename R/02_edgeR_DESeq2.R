# Multi-method sensitivity: edgeR-QLF + DESeq2 per cohort (patient-level) -> TableS3
# (Within-subtype S5 and per-source counts S5b are produced by 01_pseudobulk_meta.R.)
# Working dir: arg1, else $AD_ATLAS_DIR, else "." ; expects ec_subset.rds there.
suppressMessages({library(Seurat);library(limma);library(edgeR);library(Matrix)})
has_deseq<-requireNamespace("DESeq2",quietly=TRUE)
args<-commandArgs(trailingOnly=TRUE); setwd(if(length(args)>=1) args[[1]] else Sys.getenv("AD_ATLAS_DIR","."))
ec<-readRDS("ec_subset.rds"); md<-ec@meta.data
vec<-rownames(md)[grepl("VEC",md$celltype)]; md<-md[vec,]; ec<-subset(ec,cells=vec)
cnt<-GetAssayData(ec,slot="counts")
md$patient<-as.character(md$orig.ident)
oa<-md$Source=="Own"&md$Group=="AD"; md$patient[oa]<-sub("-[ABC]$","",md$patient[oa])
genes<-c("SELE","SELP","VCAM1","ICAM1"); cohorts<-c("GSE189795","GSE213740","GSE222318","Own")
ctrl<-c("Control"); grpof<-function(x) ifelse(x=="AD","AD","Ctrl")   # Group labels: AD, Control, AA
pbagg<-function(cells,g){f<-factor(g);M<-sparseMatrix(i=seq_along(f),j=as.integer(f),x=1,dims=c(length(f),nlevels(f)));colnames(M)<-levels(f);as.matrix(cnt[,cells,drop=FALSE]%*%M)}
S3<-data.frame()
for(co in cohorts){
 cc<-rownames(md)[md$Source==co & (md$Group=="AD"|md$Group%in%ctrl)]
 sub<-md[cc,]; grp<-grpof(sub$Group); pb<-pbagg(cc,sub$patient)
 pg<-tapply(grp,sub$patient,function(x)x[1])[colnames(pb)]; pg<-factor(pg,levels=c("Ctrl","AD"))
 if(length(unique(pg))<2||min(table(pg))<2) next
 d<-DGEList(pb); d<-calcNormFactors(d); k<-filterByExpr(d,group=pg); k[genes[genes%in%rownames(d)]]<-TRUE
 d<-d[k,,keep.lib.sizes=FALSE]; des<-model.matrix(~pg)
 tryCatch({d<-estimateDisp(d,des); tt<-topTags(glmQLFTest(glmQLFit(d,des),coef=2),n=Inf)$table
  for(g in genes) if(g%in%rownames(tt)) S3<<-rbind(S3,data.frame(cohort=co,gene=g,method="edgeR_QLF",logFC=round(tt[g,"logFC"],3),p=signif(tt[g,"PValue"],3)))},error=function(e)cat("edgeR",co,conditionMessage(e),"\n"))
 if(has_deseq){suppressMessages(library(DESeq2)); tryCatch({
  dds<-DESeqDataSetFromMatrix(round(pb),DataFrame(pg=pg),~pg); dds<-DESeq(dds,quiet=TRUE); rs<-results(dds,contrast=c("pg","AD","Ctrl"))
  for(g in genes) if(g%in%rownames(rs)) S3<<-rbind(S3,data.frame(cohort=co,gene=g,method="DESeq2",logFC=round(rs[g,"log2FoldChange"],3),p=signif(rs[g,"pvalue"],3)))},error=function(e)cat("DESeq2",co,conditionMessage(e),"\n"))}
}
write.csv(S3,"TableS3_edgeR_DESeq2.csv",row.names=FALSE); cat("DONE 02 (TableS3)\n")
