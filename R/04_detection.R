# Per-cohort per-gene detection fraction + mean logCPM (by AD/Ctrl) -> TableS7
suppressMessages({library(Seurat);library(edgeR);library(Matrix)})
args<-commandArgs(trailingOnly=TRUE); setwd(if(length(args)>=1) args[[1]] else Sys.getenv("AD_ATLAS_DIR","."))
ec<-readRDS("ec_subset.rds"); md<-ec@meta.data
vec<-rownames(md)[grepl("VEC",md$celltype)]; md<-md[vec,]; ec<-subset(ec,cells=vec)
cnt<-GetAssayData(ec,slot="counts")
md$patient<-as.character(md$orig.ident)
oa<-md$Source=="Own"&md$Group=="AD"; md$patient[oa]<-sub("-[ABC]$","",md$patient[oa])
genes<-c("SELE","SELP","VCAM1","ICAM1"); cohorts<-c("GSE189795","GSE213740","GSE222318","Own")
ctrl<-c("Control"); grpof<-function(x) ifelse(x=="AD","AD","Ctrl")   # Group labels: AD, Control, AA
pbagg<-function(cells,g){f<-factor(g);M<-sparseMatrix(i=seq_along(f),j=as.integer(f),x=1,dims=c(length(f),nlevels(f)));colnames(M)<-levels(f);as.matrix(cnt[,cells,drop=FALSE]%*%M)}
S7<-data.frame()
for(co in cohorts){
 cc<-rownames(md)[md$Source==co & (md$Group=="AD"|md$Group%in%ctrl)]
 sub<-md[cc,]; grp<-grpof(sub$Group); pb<-pbagg(cc,sub$patient)
 pg<-tapply(grp,sub$patient,function(x)x[1])[colnames(pb)]; lc<-cpm(DGEList(pb),log=TRUE,prior.count=1)
 for(gp in c("AD","Ctrl")){ cols<-which(pg==gp); if(length(cols)==0) next
  for(g in genes){ if(!g%in%rownames(pb)) next; raw<-pb[g,cols]
   S7<-rbind(S7,data.frame(cohort=co,group=gp,gene=g,n_pat=length(cols),n_nonzero=sum(raw>0),
     frac_nonzero=round(mean(raw>0),3),mean_logCPM=round(mean(lc[g,cols]),2)))}}
}
write.csv(S7,"TableS7_detection.csv",row.names=FALSE); cat("DONE 04\n")
