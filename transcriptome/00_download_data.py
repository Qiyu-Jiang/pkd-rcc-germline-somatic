"""
Download raw data for the PKD-RCC multi-layer integrative study.
  1. TCGA-KIRC expression matrix (UCSC Xena, log2(norm_count+1), 533 tumor + 72 normal)
  2. GSE289843 normalized counts (GEO supplementary) + series matrix (sample metadata)

Usage:  python 00_download_data.py
Output: ./data/kirc_xena_hiseq.tsv.gz
        ./data/GSE289843_global.normalized_counts.csv.gz
        ./data/GSE289843_series_matrix.txt.gz
"""
import os, urllib.request

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)

SOURCES = {
    # UCSC Xena: TCGA-KIRC HiSeqV2 (Illumina HiSeq, log2(x+1) normalized, gene-symbol rows)
    "kirc_xena_hiseq.tsv.gz":
        "https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/HiSeqV2.gz",
    # GEO GSE289843 (Decuypere & Borras 2025; GEO submitters): normalised counts (linear scale), supplementary file
    "GSE289843_global.normalized_counts.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289843/suppl/GSE289843_global.normalized_counts.csv.gz",
    # Series matrix: only metadata is used (Sample_title / library name mapping);
    # NOTE: the expression table inside the series matrix is EMPTY for this series.
    "GSE289843_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289843/matrix/GSE289843_series_matrix.txt.gz",
}

for fname, url in SOURCES.items():
    dest = os.path.join(DATA, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"[skip] {fname} already exists")
        continue
    print(f"[get ] {url}")
    urllib.request.urlretrieve(url, dest)
    print(f"       -> {dest} ({os.path.getsize(dest)/1e6:.1f} MB)")
print("Done.")
