# KLF2 covariate sensitivity (exploratory) -> TableS6_KLF2_sensitivity.csv (wide: base vs +KLF2)
# Per cohort/gene: AD coefficient with vs without KLF2 log-CPM covariate; plus deltaKLF2 (AD-Ctrl).
suppressMessages({library(Seurat);library(limma);library(edgeR);library(Matrix)})
args<-commandArgs(trailingOnly=TRUE); setwd(if(length(args)>=1) args[[1]] else Sys.getenv("AD_ATLAS_DIR","."))
ec<-readRDS("ec_subset.rds"); md<-ec@meta.data
vec<-rownames(md)[grepl("VEC",md$celltype)]; md<-md[vec,]; ec<-subset(ec,cells=vec)
cnt<-GetAssayData(ec,slot="counts")
md$patient<-as.character(md$orig.ident)
oa<-md$Source=="Own"&md$Group=="AD"; md$patient[oa]<-sub("-[ABC]$","",md$patient[oa])
genes<-c("SELE","SELP","VCAM1","ICAM1"); cohorts<-c("GSE189795","GSE213740","GSE222318","Own")
ctrl<-c("Control"); grpof<-function(x) ifelse(x=="AD","AD","Ctrl")   # Group labels: AD, Control, AA
pbagg<-function(cells,g){f<-factor(g);M<-sparseMatrix(i=seq_along(f),j=as.integer(f),x=1,dims=c(length(f),nlevels(f)));colnames(M)<-levels(f);as.matrix(cnt[,cells,drop=FALSE]%*%M)}
S6<-data.frame()
for(co in cohorts){
 cc<-rownames(md)[md$Source==co & (md$Group=="AD"|md$Group%in%ctrl)]
 sub<-md[cc,]; grp<-grpof(sub$Group); pb<-pbagg(cc,sub$patient)
 pg<-tapply(grp,sub$patient,function(x)x[1])[colnames(pb)]; pg<-factor(pg,levels=c("Ctrl","AD"))
 if(length(unique(pg))<2||min(table(pg))<2) next
 d<-DGEList(pb); d<-calcNormFactors(d); k<-filterByExpr(d,group=pg); k[genes[genes%in%rownames(d)]]<-TRUE
 d<-d[k,,keep.lib.sizes=FALSE]; lc<-cpm(d,log=TRUE,prior.count=1)
 klf2<-if("KLF2"%in%rownames(lc)) lc["KLF2",] else rep(0,ncol(lc))
 dk<-round(mean(klf2[pg=="AD"])-mean(klf2[pg=="Ctrl"]),2)
 db<-model.matrix(~pg); fb<-eBayes(lmFit(voom(d,db),db))
 da<-model.matrix(~pg+klf2); fa<-eBayes(lmFit(voom(d,da),da))
 for(g in genes) if(g%in%rownames(fb)) S6<-rbind(S6,data.frame(Cohort=co,nAD=sum(pg=="AD"),nCtrl=sum(pg=="Ctrl"),
   deltaKLF2_log2FC=dk,Gene=g,AD_coef_base=round(fb$coefficients[g,"pgAD"],3),AD_coef_plusKLF2=round(fa$coefficients[g,"pgAD"],3)))
}
write.csv(S6,"TableS6_KLF2_sensitivity.csv",row.names=FALSE); cat("DONE 03 (TableS6)\n")
