#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_eqtlgen_recalibration.py

Recalibrate whole-blood cis-eQTL instruments (whole-blood layer, Table 3).

The eQTLGen consortium reports Z-scores rather than effect sizes. The original
SE = 1/sqrt(N) approximation (see Supplementary_Code_MR_final.R, Section 9) is
replaced here by the standard per-SNP conversion:

    beta = Z / sqrt(2 * p * (1 - p) * N)
    SE   = 1 / sqrt(2 * p * (1 - p) * N)

where p is the effect-allele frequency in a European reference population
(dbSNP / 1000 Genomes) and N is the per-SNP sample size (eQTLGen NrSamples).

Running this script regenerates results/M3_recalibrated_IV_table.csv, the
per-SNP table underlying the recalibrated whole-blood estimates reported in
the manuscript (PKD2 IVW: OR = 1.10, P = 0.299; PKD1 Wald ratio: OR = 0.79,
P = 0.610).

Usage:
    python 04_eqtlgen_recalibration.py                 # regenerate the submitted table
    python 04_eqtlgen_recalibration.py --affile FILE   # override AF with a SNP->AF file
                                                       # (e.g., the eQTLGen SNP_AF file:
                                                       #  2018-07-18_SNP_AF_for_AlleleB_...txt.gz, cols SNP, MAF)
    python 04_eqtlgen_recalibration.py --out FILE      # custom output path

Notes:
  - rs114796149 (PKD2) is absent from the FinnGen R5 release and was therefore
    excluded from the recalibrated analysis; it is retained in the comments
    below for completeness but is not part of the output table.
  - Outcome-side values (FinnGen R5) are carried through and match the
    harmonised dataset used in the MR analysis.
"""
import argparse
import csv
import math

# Instruments as used in the analysis.
# beta_zscale / se_zscale: exposure values in the original pipeline
#   (beta = Z-score, se = 1/sqrt(N)); af_ea: effect-allele frequency
#   (dbSNP European / 1000G EUR); outcome_*: FinnGen R5 values.
ROWS = [
    # rsid,        gene,  ea, nea, beta_zscale,           se_zscale,             af_ea,  outcome_beta, outcome_se, outcome_p
    ("rs13333747", "PKD1", "C", "T", -0.0689373314904503, 0.00740314344982767, 0.16828, 0.0311, 0.0609, 0.6093),
    ("rs35382133", "PKD2", "C", "T",  0.278031942975026,  0.00576649105726268, 0.01963, 0.1604, 0.1369, 0.2411),
    ("rs75245952", "PKD2", "T", "C", -0.0745731820065074, 0.0058794038069433,  0.03317, 0.0511, 0.1189, 0.6675),
    # excluded: rs114796149 (PKD2): beta_zscale = 0.0682413477599718, se_zscale = 0.00586492611060735, af_ea = 0.01979
]

def load_af(path):
    """Load SNP -> MAF mapping from a two-column file (rsid, freq); gz supported."""
    import gzip
    opener = gzip.open if path.endswith('.gz') else open
    af = {}
    with opener(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    af[parts[0]] = float(parts[1])
                except ValueError:
                    continue
    return af

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--affile', help='optional SNP->AF file overriding the embedded frequencies')
    ap.add_argument('--out', default='M3_recalibrated_IV_table.csv')
    a = ap.parse_args()

    af_override = load_af(a.affile) if a.affile else {}
    if a.affile and not af_override:
        raise SystemExit('No frequencies could be parsed from --affile')

    out_rows = []
    for rsid, gene, ea, nea, beta_z, se_z, af_ea, ob, ose, op in ROWS:
        p = af_override.get(rsid, af_ea)
        z = beta_z / se_z                 # full-precision Z-score
        n_full = 1.0 / se_z ** 2          # per-SNP sample size (equivalent to eQTLGen NrSamples)
        var = 2.0 * p * (1.0 - p) * n_full
        beta_r = z / math.sqrt(var)
        se_r = 1.0 / math.sqrt(var)
        out_rows.append(dict(
            rsid=rsid, gene=gene, ea=ea, N=int(round(n_full)), af_ea=p, Z=round(z, 3),
            beta_recalc=round(beta_r, 6), se_recalc=round(se_r, 6),
            beta_orig_zscale=beta_z, se_orig=se_z,
            outcome_beta=ob, outcome_se=ose, outcome_p=op,
        ))

    with open(a.out, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f'Wrote {len(out_rows)} rows -> {a.out}')
    for r in out_rows:
        print(f"  {r['rsid']} ({r['gene']}): beta={r['beta_recalc']}, SE={r['se_recalc']}")
    print('Note: rs114796149 (PKD2) excluded (absent from FinnGen R5); not included in the output table.')

if __name__ == '__main__':
    main()
