"""
GSE289843 ADPKD cyst-layer differential expression: 47 cysts (4 PKD1-mutated
patients) vs 4 healthy kidney samples, normalised counts (linear scale).
Reproduces Supplementary Table S10 (primary analysis) and Table S5
(patient-level and 40-cyst sensitivity analyses).

Sample classification (ADPKD cyst analysis / GEO metadata):
  - Healthy    : sample title/source contains 'healthy'          -> n = 4
  - Cyst       : contains 'cyst' but not 'microcystic'          -> n = 47
  - Microcystic: contains 'microcystic' (excluded from primary) -> n = 4
  - The 7 cyst samples whose title ends with '2' are second sections of
    previously sampled cysts sequenced in a later batch (retained in the primary
    analysis; a sensitivity analysis restricted to the 40 individually dissected
    cysts is also provided).

Statistics (Methods):
  - Mann-Whitney U, two-sided
  - log2FC = log2(mean_cyst + 0.01) - log2(mean_healthy + 0.01)   [log2 scale of normalised counts]
  - BH-FDR over the 43-gene panel, full-precision p-values, capped at 1
  - Sensitivity 1 (S5): patient-level - per-patient means of the 47 cysts
    (4 patient-level means) vs 4 healthy
  - Sensitivity 2 (S5): 40-cyst cohort - titles not ending in '2' vs healthy

Usage:  python 02_gse289843_de.py   (requires 00_download_data.py run first)
Output: results/S10_adpkd_de.tsv
        results/S5_sensitivity.tsv
"""
import os, gzip, io, urllib.request
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data"); OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289843/"
CNT_URL = BASE + "suppl/GSE289843_global.normalized_counts.csv.gz"
SM_URL  = BASE + "matrix/GSE289843_series_matrix.txt.gz"

# HGNC current symbol -> legacy symbols used by this GEO release
# (note: this GEO matrix already uses CXCL8 natively; only the other three need mapping)
GEO_ALIAS = {"CGAS": "MB21D1", "STING1": "TMEM173", "H2AX": "H2AFX"}
KEY_GENES = [
    "PKD1","PKD2","VHL",
    "MTOR","RPTOR","RICTOR","TSC1","TSC2",
    "CTNNB1","APC","GSK3B","AXIN1",
    "HIF1A","EPAS1","VEGFA","VEGFB",
    "TNF","IL6","IL1B","CCL2","CXCL8",
    "STING1","TBK1","IRF3","CGAS",
    "H2AX","BRCA1","TP53","ATM",
    "MYC","CCND1","EGFR","MET",
    "CDKN1A","CDKN2A","PTEN",
    "PBRM1","BAP1","SETD2",
    "HAVCR1","LCN2","NPHS1","NPHS2",
]

def bh_fdr(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order]
    fdr = ranked * n / np.arange(1, n + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    out = np.empty(n); out[order] = fdr
    return np.clip(out, 0, 1)

def fetch(url, dest):
    if not os.path.exists(dest):
        print(f"[dl ] {url}")
        urllib.request.urlretrieve(url, dest)
    return dest

# --- 1. metadata from series matrix -------------------------------------------
sm = fetch(SM_URL, os.path.join(DATA, "GSE289843_series_matrix.txt.gz"))
titles = gsm = libs = None
with gzip.open(sm, "rt") as f:
    for line in f:
        line = line.rstrip("\n")
        if   line.startswith("!Sample_title"):          titles = [s.strip('"') for s in line.split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):  gsm    = [s.strip('"') for s in line.split("\t")[1:]]
        elif line.startswith("!Sample_description") and "Library name:" in line:
            libs = [s.strip('"').replace("Library name: ", "") for s in line.split("\t")[1:]]
assert titles and gsm and libs and len(titles) == len(gsm) == len(libs), "metadata parse failed"

meta = pd.DataFrame({"library": libs, "title": titles, "gsm": gsm})
def classify(title):
    t = title.lower()
    if "healthy" in t: return "Healthy"
    if "microcystic" in t: return "Microcystic"
    if "cyst" in t: return "Cyst"
    return "Other"
meta["type"] = meta["title"].map(classify)
print("[ok ] sample types:", dict(meta["type"].value_counts()))

# --- 2. expression matrix -------------------------------------------------------
cnt = fetch(CNT_URL, os.path.join(DATA, "GSE289843_global.normalized_counts.csv.gz"))
expr = pd.read_csv(cnt, index_col=0)
common = [c for c in expr.columns if c in set(meta["library"])]
expr = expr[common]; meta = meta[meta["library"].isin(common)].set_index("library").loc[common]
print(f"[ok ] matrix {expr.shape[0]} genes x {expr.shape[1]} samples")

resolve = lambda g: GEO_ALIAS.get(g, g)
missing = [g for g in KEY_GENES if resolve(g) not in expr.index]
assert not missing, f"genes missing from GEO matrix: {missing}"

cyst    = meta[meta["type"] == "Cyst"].index.tolist()
healthy = meta[meta["type"] == "Healthy"].index.tolist()
assert len(cyst) == 47 and len(healthy) == 4, f"cohort mismatch: cyst={len(cyst)}, healthy={len(healthy)}"
print(f"[ok ] cohort: {len(cyst)} cysts (4 patients) vs {len(healthy)} healthy")

# --- 3. primary DE: 47 cysts vs 4 healthy  (Table S10) --------------------------
rows = []
for g in KEY_GENES:
    a = expr.loc[resolve(g), cyst].astype(float); b = expr.loc[resolve(g), healthy].astype(float)
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    rows.append({"gene": g, "alias_in_geo": resolve(g),
                 "cyst_mean_normalized": round(float(a.mean()), 4),
                 "healthy_mean_normalized": round(float(b.mean()), 4),
                 "log2FC_cyst_vs_healthy": round(float(np.log2(a.mean() + 0.01) - np.log2(b.mean() + 0.01)), 4),
                 "n_cyst": len(a), "n_healthy": len(b),
                 "MWU_stat": int(u), "p_value": p})
s10 = pd.DataFrame(rows); s10["BH_FDR"] = bh_fdr(s10["p_value"].values)
s10 = s10.sort_values("BH_FDR").reset_index(drop=True)
s10.to_csv(os.path.join(OUT, "S10_adpkd_de.tsv"), sep="\t", index=False, float_format="%.6e")
print(f"[done] S10_adpkd_de.tsv | significant (BH<0.05): {(s10['BH_FDR']<0.05).sum()}/43")

# --- 4. sensitivity analyses (Table S5) ------------------------------------------
pat_id = meta.loc[cyst, "title"].str.extract(r"Patient\s*(\d+)")[0]
assert pat_id.nunique() == 4, f"expected 4 patients, got {pat_id.nunique()}"
pat_mean = expr[cyst].T.groupby(pat_id.values).mean().T           # 43 genes x 4 patients
p40 = [s for s in cyst if not str(meta.loc[s, "title"]).rstrip().endswith("2")]
assert len(p40) == 40, f"expected 40 first-batch cysts, got {len(p40)}"
print(f"[ok ] sensitivity: 4 patients vs 4 healthy; {len(p40)} first-batch cysts vs 4 healthy")

rows5 = []
for g in KEY_GENES:
    gene = resolve(g)
    _, p_pat = mannwhitneyu(pat_mean.loc[gene].astype(float), expr.loc[gene, healthy].astype(float), alternative="two-sided")
    _, p_40  = mannwhitneyu(expr.loc[gene, p40].astype(float),   expr.loc[gene, healthy].astype(float), alternative="two-sided")
    _, p_all = mannwhitneyu(expr.loc[gene, cyst].astype(float),  expr.loc[gene, healthy].astype(float), alternative="two-sided")
    rows5.append({"gene": g, "cyst_p": p_all, "patient_p": p_pat, "cyst40_p": p_40})
s5 = pd.DataFrame(rows5)
s5["cyst_FDR"], s5["patient_FDR"], s5["cyst40_FDR"] = bh_fdr(s5["cyst_p"]), bh_fdr(s5["patient_p"]), bh_fdr(s5["cyst40_p"])
s5 = s5[["gene","cyst_p","cyst_FDR","patient_p","patient_FDR","cyst40_p","cyst40_FDR"]] \
       .sort_values("cyst_FDR").reset_index(drop=True)
s5.to_csv(os.path.join(OUT, "S5_sensitivity.tsv"), sep="\t", index=False, float_format="%.6e")
print(f"[done] S5_sensitivity.tsv | patient sig: {(s5['patient_FDR']<0.05).sum()}/43 | "
      f"40-cyst sig: {(s5['cyst40_FDR']<0.05).sum()}/43 | concordant-direction check in S6")

print("\nSpot checks (must match the manuscript):")
chk = s10.set_index("gene")
for g, p_exp in [("PKD1", "8.00e-06"), ("MTOR", None)]:
    row = chk.loc[g]
    print(f"  {g}: cyst p={row['p_value']:.2e}, log2FC={row['log2FC_cyst_vs_healthy']}, "
          f"FDR={row['BH_FDR']:.2e}" + (f" (expect {p_exp})" if p_exp else ""))
