"""
Batch-confounding re-analysis for the GSE289843 ADPKD cyst layer.

Reconstructed batch structure (from library IDs; GC-number gap = two library
prep rounds):
  Batch 1 (GC117831-117857): 20 cysts (Patient1 x11, Patient2 x9) + 4 microcystic
  Batch 2 (GC126296-126326): 27 cysts (Patient3 x10, Patient4 x10, 7 repeat
                             dissections of Patient1/2 cysts) + ALL 4 healthy
  -> the 4 healthy controls were sequenced ONLY in Batch 2, so the primary
     47-vs-4 comparison is partially confounded with batch.

Three analyses (Supplementary Table S11):
  A. OLS with a batch covariate, all 47 cysts vs 4 healthy:
        expression ~ 1 + cyst_status + batch1 ; t-test on cyst_status
  B. Batch-2-only contrast (zero batch confounding): 27 cysts vs 4 healthy, MWU
  C. Batch effect quantification: PCA on cyst-only samples (does PC1/2 align
     with batch after removing group?)

Usage: python 03_batch_covariate_reanalysis.py
Out:   results/S11_batch_reanalysis.tsv  + console summary
"""
import os, gzip
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.api as sm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, OUT = os.path.join(HERE,"data"), os.path.join(HERE,"results")
os.makedirs(OUT, exist_ok=True)

# ---------- metadata ----------
titles = libs = chars = None
with gzip.open(os.path.join(DATA,"GSE289843_series_matrix.txt.gz"),"rt") as f:
    for line in f:
        line=line.rstrip("\n")
        if line.startswith("!Sample_title"): titles=[v.strip('"') for v in line.split("\t")[1:]]
        elif line.startswith("!Sample_description") and "Library" in line:
            libs=[v.strip('"').replace("Library name: ","") for v in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1") and "tissue:" in line:
            chars=[v.strip('"').replace("tissue: ","") for v in line.split("\t")[1:]]
meta = pd.DataFrame({"title":titles,"library":libs,"tissue":chars})
meta["gc"] = meta["library"].str.extract(r"GC(\d+)").astype(int)
meta["batch"] = np.where(meta["gc"] < 120000, "B1", "B2")
meta["group"] = np.where(meta["tissue"].str.contains("Healthy"), "healthy",
                 np.where(meta["tissue"].str.contains("microcystic"), "micro", "cyst"))
ana = meta[meta.group.isin(["cyst","healthy"])].reset_index(drop=True)   # 51 samples
print(f"[ok ] samples in model: {len(ana)}  (cyst {sum(ana.group=='cyst')}, healthy {sum(ana.group=='healthy')})")
print(f"      Batch1: {sum((ana.batch=='B1')&(ana.group=='cyst'))} cysts, 0 healthy | "
      f"Batch2: {sum((ana.batch=='B2')&(ana.group=='cyst'))} cysts + {sum((ana.batch=='B2')&(ana.group=='healthy'))} healthy")

expr = pd.read_csv(os.path.join(DATA,"GSE289843_global.normalized_counts.csv.gz"), index_col=0)
X_samples = ana.library.tolist()
KEY_GENES = ["PKD1","PKD2","VHL","MTOR","RPTOR","RICTOR","TSC1","TSC2","CTNNB1","APC","GSK3B","AXIN1",
             "HIF1A","EPAS1","VEGFA","VEGFB","TNF","IL6","IL1B","CCL2","CXCL8","STING1","TBK1","IRF3",
             "CGAS","H2AX","BRCA1","TP53","ATM","MYC","CCND1","EGFR","MET","CDKN1A","CDKN2A","PTEN",
             "PBRM1","BAP1","SETD2","HAVCR1","LCN2","NPHS1","NPHS2"]
ALIAS = {"CGAS":"MB21D1","STING1":"TMEM173","H2AX":"H2AFX"}
missing=[g for g in KEY_GENES if ALIAS.get(g,g) not in expr.index]; assert not missing, missing
Y = expr.loc[[ALIAS.get(g,g) for g in KEY_GENES], X_samples].astype(float)
Y.index = KEY_GENES

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p)
    f=p[o]*n/np.arange(1,n+1); f=np.minimum.accumulate(f[::-1])[::-1]
    return np.clip(pd.Series(f,index=o).sort_index().values,0,1)

cyst_status = (ana.group=="cyst").astype(float).values
batch1      = (ana.batch=="B1").astype(float).values
D = sm.add_constant(np.column_stack([cyst_status, batch1]))
cyst_b1  = [s for s,g,b in zip(ana.library,ana.group,ana.batch) if g=="cyst" and b=="B1"]
cyst_b2  = [s for s,g,b in zip(ana.library,ana.group,ana.batch) if g=="cyst" and b=="B2"]
healthy_s = ana[ana.group=="healthy"].library.tolist()

rows=[]
for g in KEY_GENES:
    y = Y.loc[g].values
    # A. OLS with batch covariate
    fit = sm.OLS(y, D).fit()
    beta_cyst, p_ols = fit.params[1], fit.pvalues[1]
    # log2FC-style effect from the model: adjusted mean difference (normalised-count units)
    # B. Batch2-only MWU (zero confounding)
    u_b2, p_b2 = mannwhitneyu(Y.loc[g, cyst_b2], Y.loc[g, healthy_s], alternative="two-sided")
    lfc_b2 = np.log2(Y.loc[g,cyst_b2].mean()+0.01) - np.log2(Y.loc[g,healthy_s].mean()+0.01)
    # unadjusted MWU (primary analysis, S10)
    u_un, p_un = mannwhitneyu(Y.loc[g, cyst_b1+cyst_b2], Y.loc[g, healthy_s], alternative="two-sided")
    lfc_un = np.log2(Y.loc[g,cyst_b1+cyst_b2].mean()+0.01) - np.log2(Y.loc[g,healthy_s].mean()+0.01)
    rows.append(dict(gene=g, unadj_P=p_un, unadj_log2FC=lfc_un, ols_beta=beta_cyst, ols_P=p_ols,
                     B2only_log2FC=lfc_b2, B2only_P=p_b2))
r = pd.DataFrame(rows)
r["unadj_FDR"], r["ols_FDR"], r["B2only_FDR"] = bh(r.unadj_P), bh(r.ols_P), bh(r.B2only_P)
dir_un  = np.sign(r.unadj_log2FC); dir_b2 = np.sign(r.B2only_log2FC)
r["direction_consistent"] = dir_un == dir_b2
r = r.sort_values("ols_FDR").reset_index(drop=True)
r.to_csv(os.path.join(OUT,"S11_batch_reanalysis.tsv"), sep="\t", index=False, float_format="%.6e")

# ---------- summary ----------
def nsig(col): return int((r[col]<0.05).sum())
print(f"""
===== Summary (43 genes) =====
Primary (unadjusted MWU)        significant: {nsig('unadj_FDR')}/43
A. Batch-covariate OLS          significant: {nsig('ols_FDR')}/43   | same direction: {(np.sign(r.ols_beta)==dir_un).sum()}/43
B. Batch2-only contrast         significant: {nsig('B2only_FDR')}/43 | direction consistent: {r.direction_consistent.sum()}/43

Significant unadjusted but dropped in OLS (FDR >= 0.05): {r[(r.unadj_FDR<0.05)&(r.ols_FDR>=0.05)].gene.tolist() or 'none'}
Not significant unadjusted but significant in OLS:      {r[(r.unadj_FDR>=0.05)&(r.ols_FDR<0.05)].gene.tolist() or 'none'}
Direction-flipped genes:                                {r[~r.direction_consistent].gene.tolist() or 'none'}
""")
print("Key gene details:")
key_chk = ["PKD1","IL6","CXCL8","TNF","IL1B","CCL2","STING1","CGAS","TBK1","IRF3","PBRM1","BAP1","SETD2","LCN2","HAVCR1","MTOR"]
print(r[r.gene.isin(key_chk)][["gene","unadj_log2FC","unadj_P","ols_beta","ols_P","ols_FDR","B2only_log2FC","B2only_P","direction_consistent"]].to_string(index=False))

# ---------- C. PCA on cyst-only ----------
from sklearn.decomposition import PCA
cyst_all = Y[cyst_b1+cyst_b2].T                      # 47 cysts x 43 genes
Xs = (cyst_all - cyst_all.mean()) / cyst_all.std().replace(0,1)
pc = PCA(n_components=2).fit_transform(Xs)
lab_b1 = [1 if s in cyst_b1 else 0 for s in cyst_all.index]
r_pc1 = np.corrcoef(pc[:,0], lab_b1)[0,1]
print(f"\n[C] cyst-only PCA (47x43): r(PC1, batch) = {r_pc1:.3f}  (|r| > 0.5 indicates batch clustering)")
pc_b1, pc_b2 = pc[np.array(lab_b1)==1], pc[np.array(lab_b1)==0]
from scipy.stats import ttest_ind
print(f"    PC1 mean: Batch1={pc_b1[:,0].mean():.2f} vs Batch2={pc_b2[:,0].mean():.2f}  (t-test p={ttest_ind(pc_b1[:,0],pc_b2[:,0]).pvalue:.2e})")
