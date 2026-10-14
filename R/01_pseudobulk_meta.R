# AUTHORITATIVE analysis that produced the reported S2A/S2B/S4/S5 (cleaned from ec_patient.R).
# Patient-level pseudobulk -> per-cohort limma-voom -> TableS2A
#   NB order: filterByExpr THEN calcNormFactors; SE = stdev.unscaled[,2]*sqrt(s2.post)
# Random-effects meta (REML; Hartung-Knapp) on FULL-PRECISION effects/SE via metafor -> TableS2B
# Leave-one-cohort-out (HK + conventional RE) -> TableS4
# Within-subtype 12 tests (~Source+Group, rank-checked) -> TableS5 ; per-source counts -> TableS5b
# Working dir: arg1, else $AD_ATLAS_DIR, else "." ; expects ec_subset.rds there.
suppressMessages({library(edgeR);library(limma);library(metafor);library(Seurat);library(Matrix)})
args<-commandArgs(trailingOnly=TRUE); setwd(if(length(args)>=1) args[[1]] else Sys.getenv("AD_ATLAS_DIR","."))
ec<-readRDS("ec_subset.rds")
cnt<-tryCatch(GetAssayData(ec,layer="counts"),error=function(e)GetAssayData(ec,"counts")); md<-ec@meta.data
md$patient<-as.character(md$orig.ident)
oa<-md$Source=="Own" & md$Group=="AD"; md$patient[oa]<-sub("-[ABC]$","",md$patient[oa])
adh<-c("SELE","SELP","VCAM1","ICAM1"); cohorts<-c("GSE189795","GSE213740","GSE222318","Own")
vec<-rownames(md)[md$celltype=="VEC"]; pats<-unique(md[vec,"patient"])
pb<-sapply(pats,function(p){cc<-vec[md[vec,"patient"]==p];Matrix::rowSums(cnt[,cc,drop=FALSE])}); colnames(pb)<-pats
pin<-md[match(pats,md$patient),c("patient","Source","Group")]

## ---- per-cohort limma-voom (filter THEN normalize; vectorized SE) ----
rows<-list()
for(co in cohorts){ pp<-pats[pin$Source==co & pin$Group%in%c("AD","Control")]
  grp<-factor(pin$Group[match(pp,pin$patient)],levels=c("Control","AD"))
  y<-DGEList(pb[,pp,drop=FALSE]); keep<-filterByExpr(y,group=grp); keep[adh]<-TRUE
  y<-y[keep,,keep.lib.sizes=FALSE]; y<-calcNormFactors(y); des<-model.matrix(~grp)
  v<-voom(y,des); fit<-eBayes(lmFit(v,des)); lfc<-fit$coefficients[,2]; se<-fit$stdev.unscaled[,2]*sqrt(fit$s2.post)
  for(g in adh) if(g%in%rownames(y)) rows[[paste(co,g)]]<-data.frame(Gene=g,Cohort=co,
    n_AD=sum(grp=="AD"),n_Ctrl=sum(grp=="Control"),log2FC=lfc[g],SE=se[g]) }
df<-do.call(rbind,rows)                                   # FULL precision kept in df
S2Ao<-df; S2Ao$log2FC<-round(S2Ao$log2FC,3); S2Ao$SE<-round(S2Ao$SE,3)
write.csv(S2Ao,"TableS2A_percohort_voom.csv",row.names=FALSE)   # display rounded only
S1<-unique(df[,c("Cohort","n_AD","n_Ctrl")]); write.csv(S1,"TableS1_cohort_counts.csv",row.names=FALSE)

## ---- meta + leave-one-out (metafor, full precision) ----
loo<-function(d,g,om){ hk<-rma(yi=log2FC,sei=SE,data=d,method="REML",test="knha"); zz<-rma(yi=log2FC,sei=SE,data=d,method="REML")
  data.frame(Gene=g,Omitted=om,k=nrow(d),Pooled_log2FC=round(hk$b,3),HK_CI_low=round(hk$ci.lb,3),HK_CI_high=round(hk$ci.ub,3),
    p_HartungKnapp=signif(hk$pval,4),p_conventionalRE=signif(zz$pval,4),direction=ifelse(hk$b<0,"neg","POS")) }
S2B<-data.frame(); S4<-data.frame()
for(g in adh){ d<-df[df$Gene==g,]
  hk<-rma(yi=log2FC,sei=SE,data=d,method="REML",test="knha")
  S2B<-rbind(S2B,data.frame(Gene=g,Pooled_log2FC=round(hk$b,3),KNHA_95CI_low=round(hk$ci.lb,3),
    KNHA_95CI_high=round(hk$ci.ub,3),`I2(%)`=round(hk$I2,1),KNHA_p=signif(hk$pval,4),check.names=FALSE))
  S4<-rbind(S4,loo(d,g,"(full)")); for(i in 1:nrow(d)) S4<-rbind(S4,loo(d[-i,],g,paste0("-",d$Cohort[i]))) }
S2B$BH_FDR<-signif(p.adjust(S2B$KNHA_p,"BH"),4)
write.csv(S2B,"TableS2B_meta_pooled.csv",row.names=FALSE); write.csv(S4,"TableS4_leave_one_out.csv",row.names=FALSE)

## ---- within-subtype 12 tests (~Source+Group) + rank check + per-source counts ----
S5<-data.frame(); SRC<-data.frame()
for(st in c("VEC_C0_ACKR1+POSTN+","VEC_C1_ACKR1+MYC+","VEC_C4_BTNL9+")){
  cl<-rownames(md)[md$celltype=="VEC" & md$subtype_2==st & md$Source%in%cohorts & md$Group%in%c("AD","Control")]
  p2<-unique(md[cl,"patient"]); m<-sapply(p2,function(p){cc<-cl[md[cl,"patient"]==p]; if(length(cc)<3) rep(NA_real_,nrow(cnt)) else Matrix::rowSums(cnt[,cc,drop=FALSE])})
  ok<-!is.na(m[1,]); m<-m[,ok,drop=FALSE]; p2<-p2[ok]
  src<-factor(md[match(p2,md$patient),"Source"]); grp<-factor(md[match(p2,md$patient),"Group"],levels=c("Control","AD"))
  for(s in levels(src)) SRC<-rbind(SRC,data.frame(subtype=st,Source=s,n_AD=sum(src==s&grp=="AD"),n_Ctrl=sum(src==s&grp=="Control")))
  y<-DGEList(m); keep<-filterByExpr(y,group=grp); keep[adh]<-TRUE; y<-y[keep,,keep.lib.sizes=FALSE]; y<-calcNormFactors(y)
  des_full<-if(nlevels(src)>1) model.matrix(~src+grp) else model.matrix(~grp)
  if(qr(des_full)$rank < ncol(des_full)){ des<-model.matrix(~grp); flag<-"RANK-DEFICIENT: Source dropped -> ~Group"
  } else { des<-des_full; flag<-if(nlevels(src)>1) "full-rank ~Source+Group" else "single-source ~Group" }
  cat("[",st,"]",flag," rank=",qr(des_full)$rank,"/",ncol(des_full),"\n")
  v<-voom(y,des); fit<-eBayes(lmFit(v,des)); gc<-grep("grpAD",colnames(fit$coefficients))
  for(g in adh) if(g%in%rownames(y)) S5<-rbind(S5,data.frame(subtype=st,gene=g,model=flag,
    logFC=round(fit$coefficients[g,gc],3),p=signif(fit$p.value[g,gc],3))) }
S5$BH_FDR<-signif(p.adjust(S5$p,"BH"),3)
write.csv(S5,"TableS5_subtype_12tests.csv",row.names=FALSE); write.csv(SRC,"TableS5b_subtype_source_counts.csv",row.names=FALSE)
cat("DONE 01 (S2A/S2B/S4/S5/S5b; meta on full precision; within-subtype rank-checked)\n")
