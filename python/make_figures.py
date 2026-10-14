"""Regenerate all figures from the VERIFIED R metafor result tables + fig1_coords.csv.

IMPORTANT: this script is for PLOTTING ONLY. All pooled effect sizes, confidence
intervals, I2 and p-values are read directly from the R metafor output tables
(TableS2B meta; TableS4 leave-one-out). No meta-analysis is recomputed here, so
figures are byte-consistent with the authoritative R results. (meta.py is a
separate, approximate independent cross-check and is NOT the source of any
reported value.)

Inputs (pass dir with --dir; default current):
  TableS2A_percohort_voom.csv   per-cohort log2FC/SE (squares + per-cohort CI)
  TableS2B_meta_pooled.csv      R metafor pooled log2FC/CI/I2/p   (Fig2 diamonds, I2)
  TableS3_edgeR_DESeq2.csv      multi-method (Fig4B)
  TableS4_leave_one_out.csv     R metafor leave-one-out pooled/CI/p  (Fig4A)
  TableS5_subtype_12tests.csv   within-subtype (Fig3)
  TableS6_KLF2_sensitivity.csv  KLF2 base vs +KLF2 (FigS1)
  fig1_coords.csv               UMAP coords (Fig1)
Outputs: Fig1-4 + FigS1 (pdf+png).
"""
import argparse, csv, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

GEN=["SELE","SELP","VCAM1","ICAM1"]; COH=["GSE189795","GSE213740","GSE222318","Own"]
GC={"SELE":"#c0392b","SELP":"#2980b9","VCAM1":"#27ae60","ICAM1":"#8e44ad"}
_n=lambda c:c.split("(")[0]  # 'Own(23pt)' -> 'Own'

def _f(x):
    try: return float(x)
    except: return float("nan")

def load(d):
    S2A={}
    for r in csv.DictReader(open(os.path.join(d,"TableS2A_percohort_voom.csv"))):
        S2A[(r["Gene"],_n(r["Cohort"]))]=(float(r["log2FC"]),float(r["SE"]),int(r["n_AD"]),int(r["n_Ctrl"]))
    S2B={}
    for r in csv.DictReader(open(os.path.join(d,"TableS2B_meta_pooled.csv"))):
        S2B[r["Gene"]]=dict(pooled=_f(r.get("Pooled_log2FC")),lo=_f(r.get("KNHA_95CI_low")),
            hi=_f(r.get("KNHA_95CI_high")),I2=_f(r.get("I2(%)",r.get("I2"))),p=_f(r.get("KNHA_p")),fdr=_f(r.get("BH_FDR")))
    S3={}
    for r in csv.DictReader(open(os.path.join(d,"TableS3_edgeR_DESeq2.csv"))):
        S3[(_n(r["cohort"]),r["gene"],r["method"])]=float(r["logFC"])
    S4={}
    for r in csv.DictReader(open(os.path.join(d,"TableS4_leave_one_out.csv"))):
        g=r["Gene"]; om=r.get("Omitted",r.get("Omitted_cohort"))
        S4[(g,om)]=dict(pooled=_f(r.get("Pooled_log2FC",r.get("Pooled"))),
            lo=_f(r.get("HK_CI_low",r.get("CI_low"))),hi=_f(r.get("HK_CI_high",r.get("CI_high"))),
            p=_f(r.get("p_HartungKnapp",r.get("KNHA_p",r.get("p_HK")))))
    S5={}
    for r in csv.DictReader(open(os.path.join(d,"TableS5_subtype_12tests.csv"))):
        S5[(r["subtype"],r["gene"])]=(float(r["logFC"]),float(r["BH_FDR"]))
    return S2A,S2B,S3,S4,S5

def fig2(S2A,S2B,out):
    fig=plt.figure(figsize=(11,10)); gs=GridSpec(1,3,width_ratios=[2.6,5.4,2.3],wspace=0.02)
    axL,axM,axR=[fig.add_subplot(gs[i]) for i in range(3)]
    for ax in (axL,axR): ax.set_xlim(0,1); ax.axis("off")
    axM.set_xlim(-6,0.7)
    rows=[]
    for g in GEN:
        rows.append(("hdr",g)); [rows.append(("coh",(g,c))) for c in COH]; rows.append(("pool",g)); rows.append(("gap",None))
    rows=rows[:-1]; n=len(rows); yp={i:n-i for i in range(n)}
    ymin=min(yp.values())-1; ymax=max(yp.values())+1
    for ax in (axL,axM,axR): ax.set_ylim(ymin,ymax)
    axM.axvline(0,color="k",ls="--",lw=.9)
    for i,(t,dat) in enumerate(rows):
        y=yp[i]
        if t=="hdr":
            b=S2B[dat]; axL.text(0,y,dat,fontsize=11,fontweight="bold",va="center")
            axL.text(0,y-0.55,f"I²={b['I2']:.0f}%  KNHA p={b['p']:.3f}  FDR={b['fdr']:.3f}",fontsize=7.5,color="#444",va="center")
        elif t=="coh":
            g,c=dat; lfc,se,na,nc=S2A[(g,c)]; lo=lfc-1.96*se; hi=lfc+1.96*se
            axL.text(0,y,c,fontsize=9,va="center"); axL.text(1,y,f"AD={na}/Ctrl={nc}",fontsize=8,color="#555",ha="right",va="center")
            axM.plot([lo,hi],[y,y],color="#555",lw=1.5); axM.scatter(lfc,y,s=55,marker="s",color="#34495e",zorder=3)
            axR.text(0,y,f"{lfc:.2f} [{lo:.2f}, {hi:.2f}]",fontsize=8.5,va="center")
        elif t=="pool":
            b=S2B[dat]; mu,lo,hi=b["pooled"],b["lo"],b["hi"]; h=.34
            axL.text(0,y,"Pooled (RE, KNHA)",fontsize=9.5,fontweight="bold",color="#b4231f",va="center")
            axM.fill([lo,mu,hi,mu],[y,y+h,y,y-h],color="#b4231f",zorder=3)
            axR.text(0,y,f"{mu:.2f} [{lo:.2f}, {hi:.2f}]",fontsize=9,fontweight="bold",color="#b4231f",va="center")
    axM.set_xlabel("log2 fold change (AD vs Control), patient-level",fontsize=9.5); axM.set_yticks([])
    for s in ["top","right","left"]: axM.spines[s].set_visible(False)
    fig.suptitle("Fig 2  Patient-level random-effects meta-analysis (values from R metafor)",fontsize=12,fontweight="bold",y=0.98)
    plt.subplots_adjust(top=0.93,bottom=0.06,left=0.01,right=0.99)
    for e in ("pdf","png"): fig.savefig(os.path.join(out,f"Fig2_meta_forest.{e}"),dpi=150)

def fig3(S5,out):
    subs=["VEC_C0_ACKR1+POSTN+","VEC_C1_ACKR1+MYC+","VEC_C4_BTNL9+"]; lab=["C0 ACKR1+POSTN+","C1 ACKR1+MYC+","C4 BTNL9+"]
    M=np.zeros((4,3)); F=np.zeros((4,3))
    for i,g in enumerate(GEN):
        for j,s in enumerate(subs): M[i,j],F[i,j]=S5[(s,g)]
    star=lambda p:"***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else ""
    fig,ax=plt.subplots(figsize=(5.2,5)); im=ax.imshow(M,cmap="Blues_r",vmin=-2.1,vmax=0,aspect="auto")
    for i in range(4):
        for j in range(3):
            c="white" if M[i,j]<-1.2 else "black"
            ax.text(j,i-.12,f"{M[i,j]:.2f}",ha="center",va="center",fontsize=10,color=c)
            ax.text(j,i+.22,star(F[i,j]),ha="center",va="center",fontsize=11,fontweight="bold",color=c)
    ax.set_xticks(range(3)); ax.set_xticklabels(lab,fontsize=9); ax.set_yticks(range(4)); ax.set_yticklabels(GEN,fontsize=11)
    ax.set_title("Within-subtype AD vs control (patient-level)\nlog2FC; *** <0.001 ** <0.01 * <0.05",fontsize=9.5)
    fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04); plt.tight_layout()
    for e in ("pdf","png"): fig.savefig(os.path.join(out,f"Fig3_subtype_heatmap.{e}"),dpi=150)

def fig4(S2B,S3,S4,out):
    fig=plt.figure(figsize=(12,6.6)); axB=fig.add_subplot(1,2,1); axC=fig.add_subplot(1,2,2)
    yp=0; yt=[]; yl=[]
    def getrow(g,om):  # tolerant key lookup
        for k in [(g,om),(g,om.replace("-","- ")),(g,"(full)"),(g,"(full/full4)"),(g,"(none/full4)")]:
            if k in S4: return S4[k]
        return None
    for g in GEN:
        labs=["(full)"]+[f"-{c}" for c in COH]
        for lb in labs[::-1]:
            r=S4.get((g,lb)) or (S4.get((g,"(full)")) if lb=="(full)" else None)
            if r is None:  # match alt 'full' spellings
                r=next((v for (gg,oo),v in S4.items() if gg==g and (oo==lb or (lb=="(full)" and "full" in oo))),None)
            if r is None: continue
            mk="D" if "full" in lb else "o"
            axB.plot([r["lo"],r["hi"]],[yp,yp],color=GC[g],lw=1.7); axB.scatter(r["pooled"],yp,color=GC[g],marker=mk,s=42,zorder=3)
            yl.append(f"{g} {lb}"+("  *" if (r["p"]==r["p"] and r["p"]<0.05) else "")); yt.append(yp); yp+=1
        yp+=0.7
    axB.axvline(0,color="k",lw=.8); axB.set_yticks(yt); axB.set_yticklabels(yl,fontsize=7.5)
    axB.set_xlabel("pooled log2FC [HK 95% CI] (R metafor)",fontsize=9); axB.set_title("A  Leave-one-cohort-out",fontsize=10.5,fontweight="bold",loc="left")
    for g in GEN:
        for c in COH:
            x=S2B[g]["pooled"] if False else None
    # Fig4B: per-cohort voom vs edgeR/DESeq2 — needs per-cohort voom; read from S2A via caller
    return fig,axB,axC

def fig4_full(S2A,S2B,S3,S4,out):
    fig,axB,axC=fig4(S2B,S3,S4,out)
    for g in GEN:
        for c in COH:
            x=S2A[(g,c)][0]
            if (c,g,"edgeR_QLF") in S3: axC.scatter(x,S3[(c,g,"edgeR_QLF")],color=GC[g],marker="o",s=42,alpha=.85)
            if (c,g,"DESeq2") in S3: axC.scatter(x,S3[(c,g,"DESeq2")],color=GC[g],marker="^",s=42,alpha=.85)
    lim=[-4,.3]; axC.plot(lim,lim,"k--",lw=.8,alpha=.6); axC.set_xlim(lim); axC.set_ylim(lim)
    axC.set_xlabel("limma-voom log2FC",fontsize=9); axC.set_ylabel("edgeR (o) / DESeq2 (^) log2FC",fontsize=9)
    axC.set_title("B  Multi-method concordance",fontsize=10.5,fontweight="bold",loc="left")
    axC.legend(handles=[Line2D([0],[0],marker="s",color="w",markerfacecolor=GC[g],label=g,markersize=8) for g in GEN],fontsize=8,loc="lower right",frameon=False)
    fig.suptitle("Fig 4  Robustness (R metafor leave-one-out + multi-method)",fontsize=12,fontweight="bold",y=0.99)
    plt.subplots_adjust(top=0.86,bottom=0.1,left=0.13,right=0.97,wspace=0.3)
    for e in ("pdf","png"): fig.savefig(os.path.join(out,f"Fig4_robustness.{e}"),dpi=150)

def figS1(d,out):
    fp=os.path.join(d,"TableS6_KLF2_sensitivity.csv")
    if not os.path.exists(fp): print("TableS6 absent; skipping FigS1"); return
    base={}; adj={}
    for r in csv.DictReader(open(fp)):
        c=_n(r["Cohort"]); base[(c,r["Gene"])]=float(r["AD_coef_base"]); adj[(c,r["Gene"])]=float(r["AD_coef_plusKLF2"])
    fig=plt.figure(figsize=(11,3.4))
    for gi,g in enumerate(GEN):
        ax=fig.add_subplot(1,4,gi+1)
        for i,c in enumerate(COH):
            if (c,g) not in base: continue
            b,a2=base[(c,g)],adj[(c,g)]; ax.plot([i,i],[b,a2],color="grey",lw=1,zorder=1)
            ax.scatter(i,b,facecolors="white",edgecolors=GC[g],s=50,zorder=2); ax.scatter(i,a2,color=GC[g],s=50,zorder=3)
        ax.axhline(0,color="k",lw=.6); ax.set_xticks(range(4)); ax.set_xticklabels(["95","740","2318","Own"],fontsize=7.5)
        ax.set_title(g,fontsize=10,color=GC[g],fontweight="bold")
        if gi==0: ax.set_ylabel("AD coef (log2FC)\nopen=base, filled=+KLF2",fontsize=8.5)
        ax.set_ylim(-4.2,0.5); ax.tick_params(labelsize=7.5)
    fig.suptitle("Fig S1  KLF2 covariate sensitivity (exploratory): open=base, filled=+KLF2",fontsize=10,y=1.02)
    plt.subplots_adjust(top=0.72,bottom=0.12,left=0.09,right=0.97,wspace=0.35)
    for e in ("pdf","png"): fig.savefig(os.path.join(out,f"FigS1_KLF2_sensitivity.{e}"),dpi=150,bbox_inches="tight")

def fig1(d,out):
    import pandas as pd
    fp=os.path.join(d,"fig1_coords.csv")
    if not os.path.exists(fp): print("fig1_coords.csv absent; skipping Fig1"); return
    df=pd.read_csv(fp)
    order=["VEC_C0_ACKR1+POSTN+","VEC_C1_ACKR1+MYC+","VEC_C2_ACKR1+IL6+","VEC_C3_IGFBP5+","VEC_C4_BTNL9+",
    "VEC_C5_SEMA3G+","VEC_C6_ELN+","VEC_C7_TBX3+","VEC_C8_SPP1+","VEC_C9_RGS5+","VEC_C10_CCL5+","LEC_C0_CCL21+","LEC_C1_ANGPT2+"]
    order=[s for s in order if s in set(df['subtype'])]
    pal=plt.get_cmap("tab20").colors; sc={s:pal[i%20] for i,s in enumerate(order)}
    fig=plt.figure(figsize=(15,6.4)); gs=GridSpec(2,3,height_ratios=[6,1.3],hspace=.05,wspace=.12)
    def um(ax,t):
        ax.set_xticks([]);ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.set_xlabel("UMAP1",fontsize=9); ax.set_title(t,fontsize=12,fontweight="bold",loc="left")
    axA=fig.add_subplot(gs[0,0])
    for s in order: dd=df[df['subtype']==s]; axA.scatter(dd.UMAP1,dd.UMAP2,s=2,c=[sc[s]],linewidths=0,rasterized=True)
    um(axA,"A  Endothelial subtypes (13)"); axA.set_ylabel("UMAP2",fontsize=9)
    axB=fig.add_subplot(gs[0,1])
    for k,c in {"VEC":"#2457d6","LEC":"#d6294a"}.items():
        dd=df[df['celltype']==k]; axB.scatter(dd.UMAP1,dd.UMAP2,s=2,c=c,linewidths=0,rasterized=True,label=f"{k} ({len(dd):,} shown)")
    um(axB,"B  VEC vs LEC"); axB.legend(markerscale=4,fontsize=9,loc="upper right",frameon=False)
    axC=fig.add_subplot(gs[0,2])
    for k,c in {"Control":"#9aa0a6","AA":"#e8a33d","AD":"#b4231f"}.items():
        dd=df[df['Group']==k]; axC.scatter(dd.UMAP1,dd.UMAP2,s=2,c=c,linewidths=0,rasterized=True,label=k)
    um(axC,"C  Disease group"); axC.legend(markerscale=4,fontsize=9,loc="upper right",frameon=False)
    axL=fig.add_subplot(gs[1,:]); axL.axis("off")
    h=[Line2D([0],[0],marker='o',color='w',markerfacecolor=sc[s],markersize=8,label=s) for s in order]
    axL.legend(handles=h,ncol=5,fontsize=8.5,loc="center",frameon=False,title="Endothelial subtypes (panel A)",title_fontsize=9)
    fig.suptitle("Fig 1  Endothelial cells of the integrated atlas (43,831 EC; 30,000 subsampled for display)",fontsize=12.5,fontweight="bold",y=.99)
    plt.subplots_adjust(top=.9,bottom=.02,left=.05,right=.98)
    for e in ("pdf","png"): fig.savefig(os.path.join(out,f"Fig1_endothelial_UMAP.{e}"),dpi=150)

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--dir",default="."); ap.add_argument("--out",default=".")
    a=ap.parse_args(); os.makedirs(a.out,exist_ok=True)
    S2A,S2B,S3,S4,S5=load(a.dir)
    fig2(S2A,S2B,a.out); fig3(S5,a.out); fig4_full(S2A,S2B,S3,S4,a.out); figS1(a.dir,a.out); fig1(a.dir,a.out)
    print("figures written to",a.out,"(all meta values read from R metafor tables S2B/S4)")
