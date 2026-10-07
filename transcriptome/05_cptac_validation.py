"""
05_cptac_validation.py

Independent validation of the manuscript's transcriptomic findings in the
CPTAC clear-cell renal cell carcinoma (CCRCC) proteogenomic cohort
(no sample overlap with TCGA-KIRC).

Two analyses (Results, "Validation in an independent cohort" section):
  A. Direction-of-effect validation for the 43-gene panel:
     protein (gene-level MS abundance) and RNA (RNA-seq log2 FPKM),
     tumour vs adjacent normal; Mann-Whitney U + BH-FDR within the panel.
  B. Mutation-expression correspondence for BAP1 / SETD2 / PBRM1 / VHL / KDM5C:
     own-gene protein and RNA abundance in mutation-positive vs wild-type tumours.

Inputs (public, downloaded automatically where possible):
  - Protein tumour matrix : HS_CPTAC_CCRCC_proteome_Tumor.cct          (LinkedOmics)
  - Protein normal matrix : HS_CPTAC_CCRCC_proteome_Normal.cct
  - RNA tumour matrix     : HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct
  - RNA normal matrix     : HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct
  - Clinical metadata     : CCRCC_meta.txt (mutation calls BIN 0/1)

Usage:  python 05_cptac_validation.py [--data-dir DIR] [--out-dir DIR]
Outputs (written to --out-dir, default ./results):
  S13_cptac_validation.tsv          # 43-gene panel, protein & RNA log2FC / P / FDR
  S13b_mutation_expression.tsv      # mutation vs wild-type own-gene abundance
"""
import argparse, os, urllib.request
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

PANEL = ["PKD1","PKD2","VHL","MTOR","RPTOR","RICTOR","TSC1","TSC2",
         "CTNNB1","APC","GSK3B","AXIN1","HIF1A","EPAS1","VEGFA","VEGFB",
         "TNF","IL6","IL1B","CCL2","CXCL8","STING1","TBK1","IRF3",
         "CGAS","H2AX","BRCA1","TP53","ATM","MYC","CCND1","EGFR","MET",
         "CDKN1A","CDKN2A","PTEN","PBRM1","BAP1","SETD2",
         "HAVCR1","LCN2","NPHS1","NPHS2"]

FILES = {
    "prot_T": "HS_CPTAC_CCRCC_proteome_Tumor.cct",
    "prot_N": "HS_CPTAC_CCRCC_proteome_Normal.cct",
    "rna_T":  "HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct",
    "rna_N":  "HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct",
    "meta":   "CCRCC_meta.txt",
}
BASE = "https://www.linkedomics.org/data_download/CPTAC-CCRCC/"
META_BASE = "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/"

def bh(p):
    p = np.asarray(p, float)
    mask = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pm = p[mask]; n = len(pm)
    o = np.argsort(pm); f = pm[o] * n / np.arange(1, n + 1)
    f = np.minimum.accumulate(f[::-1])[::-1]
    q[mask] = np.clip(f[np.argsort(o)], 0, 1)
    return q

def load(f):
    df = pd.read_csv(f, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = [str(c).strip() for c in df.columns]
    return df.apply(pd.to_numeric, errors="coerce")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="./data")
    ap.add_argument("--out-dir", default="./results")
    a = ap.parse_args()
    os.makedirs(a.data_dir, exist_ok=True); os.makedirs(a.out_dir, exist_ok=True)

    paths = {}
    for k, fname in FILES.items():
        dest = os.path.join(a.data_dir, fname)
        if not (os.path.exists(dest) and os.path.getsize(dest) > 1000):
            url = META_BASE + fname if k == "meta" else BASE + fname
            print(f"[dl ] {url}")
            urllib.request.urlretrieve(url, dest)
        paths[k] = dest

    pT, pN = load(paths["prot_T"]), load(paths["prot_N"])
    rT, rN = load(paths["rna_T"]), load(paths["rna_N"])

    # ---- A. panel direction-of-effect validation ----
    rows = []
    for g in PANEL:
        rec = {"gene": g}
        for layer, (T, N) in {"protein": (pT, pN), "RNA": (rT, rN)}.items():
            if g in T.index and g in N.index:
                t = T.loc[g].dropna(); n = N.loc[g].dropna()
                rec[f"{layer}_nT"], rec[f"{layer}_nN"] = len(t), len(n)
                rec[f"{layer}_log2FC"] = float(t.mean() - n.mean())
                rec[f"{layer}_P"] = float(mannwhitneyu(t, n, alternative="two-sided")[1])
            else:
                for c in (f"{layer}_nT", f"{layer}_nN", f"{layer}_log2FC", f"{layer}_P"):
                    rec[c] = np.nan
        rows.append(rec)
    s13 = pd.DataFrame(rows)
    s13["protein_FDR"] = bh(s13["protein_P"].values)
    s13["RNA_FDR"] = bh(s13["RNA_P"].values)
    s13.to_csv(os.path.join(a.out_dir, "S13_cptac_validation.tsv"), sep="\t", index=False, float_format="%.6e")
    print(f"[done] protein FDR<0.05: {(s13['protein_FDR']<0.05).sum()}/43 | RNA FDR<0.05: {(s13['RNA_FDR']<0.05).sum()}/43")

    # ---- B. mutation-expression correspondence ----
    meta = pd.read_csv(paths["meta"], sep="\t")
    meta = meta[meta["case_id"] != "data_type"].set_index("case_id")
    meta.index = meta.index.astype(str).str.strip()
    mrows = []
    for g in ["BAP1", "SETD2", "PBRM1", "VHL", "KDM5C"]:
        mutcol = f"{g}_mutation"
        mut = set(meta.index[pd.to_numeric(meta[mutcol], errors="coerce") == 1])
        wt  = set(meta.index[pd.to_numeric(meta[mutcol], errors="coerce") == 0])
        rec = {"gene": g, "n_mut": len(mut), "n_wt": len(wt)}
        for layer, df in (("protein", pT), ("RNA", rT)):
            if g not in df.index:
                for c in (f"{layer}_mut", f"{layer}_wt", f"{layer}_diff", f"{layer}_P"):
                    rec[c] = np.nan
                continue
            mu = df.loc[g, [c for c in df.columns if str(c).strip() in mut]].dropna()
            wt_ = df.loc[g, [c for c in df.columns if str(c).strip() in wt]].dropna()
            rec[f"{layer}_mut"], rec[f"{layer}_wt"] = float(mu.mean()), float(wt_.mean())
            rec[f"{layer}_diff"] = float(mu.mean() - wt_.mean())
            rec[f"{layer}_P"] = float(mannwhitneyu(mu, wt_, alternative="two-sided")[1])
        mrows.append(rec)
    s13b = pd.DataFrame(mrows)
    s13b["protein_FDR"] = bh(s13b["protein_P"].values)
    s13b["RNA_FDR"] = bh(s13b["RNA_P"].values)
    s13b.to_csv(os.path.join(a.out_dir, "S13b_mutation_expression.tsv"), sep="\t", index=False, float_format="%.6e")
    print("[done] S13b written")
    print(s13b[["gene","n_mut","n_wt","protein_diff","protein_P","RNA_diff","RNA_P"]].to_string(index=False))

if __name__ == "__main__":
    main()
