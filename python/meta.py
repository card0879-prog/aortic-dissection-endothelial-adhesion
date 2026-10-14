"""INDEPENDENT, APPROXIMATE cross-check of the meta-analysis (NOT the authoritative source).

The authoritative pooled effects, CIs, I2 and p-values are produced by R metafor
(R/01_pseudobulk_meta.R -> TableS2B, TableS4) and are what the manuscript and figures report.
This script re-derives the pooled effect, 95% CI and Hartung-Knapp / conventional p independently,
as a sanity check; the pooled effect/CI/p match metafor to 3 decimals. NOTE: the I2 printed here
uses a Q-statistic estimator and WILL DIFFER from metafor's REML-based I2 (e.g. SELP ~0% here vs
15.6% in metafor); use the metafor I2 (TableS2B), not this one.
Usage: python meta.py TableS2A_percohort_voom.csv"""
import sys, csv, numpy as np
from scipy import stats
def reml_hk(y,se):
    y=np.array(y,float); v=np.array(se,float)**2; k=len(y); t2=0.0
    for _ in range(300):
        w=1/(v+t2); mu=np.sum(w*y)/np.sum(w)
        new=max(0.0,np.sum(w**2*((y-mu)**2-v))/np.sum(w**2)+1/np.sum(w))
        if abs(new-t2)<1e-11: t2=new; break
        t2=new
    w=1/(v+t2); mu=np.sum(w*y)/np.sum(w); q=np.sum(w*(y-mu)**2)/(k-1); se_hk=np.sqrt(q/np.sum(w))
    tcr=stats.t.ppf(0.975,k-1); p_hk=2*stats.t.sf(abs(mu/se_hk),k-1)
    p_z=2*stats.norm.sf(abs(mu/np.sqrt(1/np.sum(w))))
    wf=1/v; muf=np.sum(wf*y)/np.sum(wf); Q=np.sum(wf*(y-muf)**2); I2=max(0,(Q-(k-1))/Q)*100
    return float(mu),float(mu-tcr*se_hk),float(mu+tcr*se_hk),float(I2),float(p_hk),float(p_z)
def main(path):
    rows=list(csv.DictReader(open(path))); genes=["SELE","SELP","VCAM1","ICAM1"]
    print("gene pooled [CI] I2_approx(NON-authoritative) p_HK  (I2 -> use metafor TableS2B)")
    for g in genes:
        d=[r for r in rows if r["Gene"]==g]
        y=[float(r["log2FC"]) for r in d]; se=[float(r["SE"]) for r in d]
        mu,lo,hi,I2,phk,pz=reml_hk(y,se)
        print(f"{g} {mu:.3f} [{lo:.3f},{hi:.3f}] I2_approx={I2:.1f} p_HK={phk:.4f}")
        for i,r in enumerate(d):  # leave-one-out
            idx=[j for j in range(len(d)) if j!=i]
            m=reml_hk([y[j] for j in idx],[se[j] for j in idx])
            print(f"   -{r['Cohort']}: {m[0]:.3f} p_HK={m[4]:.4f} p_conv={m[5]:.4f}")
if __name__=="__main__": main(sys.argv[1] if len(sys.argv)>1 else "TableS2A_percohort_voom.csv")
