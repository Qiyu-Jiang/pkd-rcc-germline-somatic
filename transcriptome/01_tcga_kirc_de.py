"""
TCGA-KIRC differential expression: 533 tumours (-01) vs 72 adjacent normals (-11).
Reproduces Supplementary Table S7 (column structure identical).

Sample classification (TCGA-KIRC analysis):
  - Tumour  : TCGA barcode field 4 starts with '01'  -> n = 533
  - Normal  : TCGA barcode field 4 starts with '11'  -> n = 72
  - Strictly primary tumours only. De novo primary tumours (code '05',
    e.g. TCGA-DV-A4W0-05) are EXCLUDED; the script asserts this.

Statistics (Methods):
  - Mann-Whitney U, two-sided
  - log2FC = mean(tumour) - mean(normal) on the Xena log2(norm_count+1) scale
  - BH-FDR over the 43-gene panel, computed on full-precision p-values,
    capped at 1.0. p-values are written with 6-digit scientific notation.

Usage:  python 01_tcga_kirc_de.py   (requires 00_download_data.py run first)
Output: results/S7_tcga_kirc_de.tsv
"""
import os, gzip, urllib.request
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data"); OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
XENA_URL = "https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/HiSeqV2.gz"

# --- 43-gene panel, HGNC current symbols -------------------------------------
# Xena matrix still uses some legacy symbols; map them explicitly.
XENA_ALIAS = {"CGAS": "C6orf150", "STING1": "TMEM173", "H2AX": "H2AFX", "CXCL8": "IL8"}
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
    "PBRM1","BAP1","SETD2",          # ccRCC epigenetic suppressor triad (TCGA 2013)
    "HAVCR1","LCN2","NPHS1","NPHS2", # tubular injury / podocyte markers
]

def bh_fdr(p):
    """Benjamini-Hochberg FDR with monotonicity enforcement and cap at 1."""
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order]
    fdr = ranked * n / np.arange(1, n + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    out = np.empty(n); out[order] = fdr
    return np.clip(out, 0, 1)

# --- load matrix --------------------------------------------------------------
path = os.path.join(DATA, "kirc_xena_hiseq.tsv.gz")
if not os.path.exists(path):
    print("[dl ] fetching KIRC matrix from Xena ...")
    urllib.request.urlretrieve(XENA_URL, path)
mat = pd.read_csv(path, sep="\t", index_col=0)

# resolve panel genes against the matrix via aliases
present, missing = [], []
for g in KEY_GENES:
    xg = XENA_ALIAS.get(g, g)
    (present if xg in mat.index else missing).append(g)
assert len(present) == 43, f"panel incomplete: {len(present)}/43, missing={missing}"
print(f"[ok ] 43/43 panel genes resolved (missing: {missing})")

def barcode_code(sid): return str(sid).split("-")[3][:2]
codes = pd.Series({s: barcode_code(s) for s in mat.columns})
tumour = codes[codes == "01"].index.tolist()
normal = codes[codes == "11"].index.tolist()

# hard guards: exact paper cohort
assert len(tumour) == 533, f"expected 533 primary tumours (-01), got {len(tumour)}"
assert len(normal) == 72,  f"expected 72 normals (-11), got {len(normal)}"
assert not any(s.startswith("TCGA-DV-A4W0-05") for s in tumour), \
    "de novo tumour TCGA-DV-A4W0-05 must be excluded (-01 only cohort)"
print(f"[ok ] cohort: {len(tumour)} tumours (-01) vs {len(normal)} normals (-11); -05 excluded")

# --- DE ------------------------------------------------------------------------
rows = []
for g in KEY_GENES:
    xg = XENA_ALIAS.get(g, g)
    a = mat.loc[xg, tumour].astype(float).dropna()
    b = mat.loc[xg, normal].astype(float).dropna()
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    rows.append({
        "gene": g,
        "alias_in_xena": xg,
        "log2FC_tumor_vs_normal": round(float(a.mean() - b.mean()), 6),
        "tumor_mean": round(float(a.mean()), 6),
        "normal_mean": round(float(b.mean()), 6),
        "tumor_n": len(a), "normal_n": len(b),
        "U_statistic": int(u),
        "p_value": p,                        # full precision internally
    })
de = pd.DataFrame(rows)
de["BH_FDR"] = bh_fdr(de["p_value"].values)

# output: 6-digit scientific notation for p / FDR (full precision retained)
for col in ("p_value", "BH_FDR"):
    de[col + "_display"] = de[col].map(lambda x: f"{x:.6e}")
de = de.sort_values("BH_FDR").reset_index(drop=True)
de.to_csv(os.path.join(OUT, "S7_tcga_kirc_de.tsv"), sep="\t", index=False,
          float_format="%.6e")

sig = (de["BH_FDR"] < 0.05).sum()
print(f"[done] S7_tcga_kirc_de.tsv  |  significant (BH<0.05): {sig}/43")
print(de[["gene","log2FC_tumor_vs_normal","U_statistic","p_value_display","BH_FDR_display"]]
      .head(10).to_string(index=False))
print("\nSpot checks (must match the manuscript):")
chk = de.set_index("gene")
for g, u_exp, p_exp in [("PKD1", 33049, "2.36e-23"), ("MTOR", 929, "2.68e-39")]:
    ok_u = chk.loc[g, "U_statistic"] == u_exp
    ok_p = f"{chk.loc[g,'p_value']:.2e}" == p_exp
    print(f"  {g}: U={chk.loc[g,'U_statistic']} ({'OK' if ok_u else 'FAIL'}), "
          f"p={chk.loc[g,'p_value']:.2e} ({'OK' if ok_p else 'FAIL'})")
