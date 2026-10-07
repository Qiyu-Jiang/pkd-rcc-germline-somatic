# PKD–RCC multi-layer integrative study: analysis code

Code and per-analysis outputs for the manuscript *"Germline exclusion and somatic convergence in ADPKD-associated renal cell carcinoma"*.

## Repository layout

```
├── Supplementary_Code_MR_final.R      # Two-sample MR: clinical layer (21 IVs) + positive control
├── rerun_21iv_local.R                 # Local re-run of the clinical-layer MR (TwoSampleMR)
├── PKD_Clinical_22_IVs_Details.csv    # Exposure instruments, harmonisation-stage snapshot (22 SNPs; see Methods notes)
├── outcome_dat_22IV_combined.csv      # Harmonised outcome data (same snapshot; 3 rows are API proxies)
├── RERUN_het_21IV.csv                 # Heterogeneity statistics (Cochran's Q)
├── RERUN_LOO_21IV.csv                 # Leave-one-out results
├── RERUN_pleiotropy_21IV.csv          # MR-Egger intercept / MR-PRESSO
├── transcriptome/
│   ├── 00_download_data.py            # Download TCGA-KIRC (Xena) and GSE289843 (GEO)
│   ├── 01_tcga_kirc_de.py             # TCGA-KIRC DE, 43-gene panel (MWU + BH)
│   ├── 02_gse289843_de.py             # GSE289843 ADPKD cyst DE (cyst vs healthy; patient-level sensitivity)
│   ├── 03_batch_covariate_reanalysis.py  # Batch-covariate sensitivity (RINT-OLS + batch-2-only + PCA)
│   ├── 04_eqtlgen_recalibration.py    # Per-SNP recalibration of whole-blood eQTL instruments (Table 3)
│   └── 05_cptac_validation.py         # Independent CPTAC-CCRCC validation (panel direction + mutation-expression; Table S13)
└── results/
    ├── Table_S11_batch_reanalysis.xlsx   # Batch-covariate sensitivity (24/43 robust (BH FDR<0.05; 25/43 at P<0.05), PCA r=0.216)
    ├── M1_covariate_expression_v2.csv    # 533 TCGA-KIRC tumours x 43 genes + stage/grade covariates
    ├── M1_covariate_adjusted_results.csv # Stage/grade-adjusted DE (beta, SE, t, P, BH q)
    ├── M3_recalibrated_IV_table.csv      # Whole-blood eQTL IVs: per-SNP recalibrated beta/SE
    │   (script: transcriptome/04_eqtlgen_recalibration.py)
    ├── S13_cptac_validation.tsv          # CPTAC-CCRCC panel validation (protein & RNA log2FC, MWU P, BH FDR)
    └── S13b_mutation_expression.tsv      # Mutation vs wild-type own-gene abundance (BAP1/SETD2/PBRM1/VHL/KDM5C)
        (script: transcriptome/05_cptac_validation.py)
```

## Environment

- R >= 4.2: TwoSampleMR, ieugwasr, MRPRESSO, data.table
- Python >= 3.9: pandas, numpy, scipy, statsmodels, openpyxl

## Data sources

| Dataset | Source | Accession / ID |
|---|---|---|
| Clinical PKD GWAS (exposure) | IEU OpenGWAS | ebi-a-GCST90018904 |
| Renal cancer outcome | FinnGen R5 | finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC (971 cases / 174,006 controls) |
| TCGA-KIRC expression | UCSC Xena | TCGA.KIRC.sampleMap / GDC (533 tumours / 72 normal) |
| ADPKD cyst RNA-seq | GEO | GSE289843 |
| Whole-blood cis-eQTL | eQTLGen consortium | cis-eQTLs full 2018-09-05 |

## Methods notes

- The whole-blood eQTL layer reports Z-scores in the source data; per-SNP exposure beta and SE were
  reconstructed as beta = Z / sqrt(2*p*(1-p)*N) and SE = 1 / sqrt(2*p*(1-p)*N), using European allele
  frequencies (dbSNP) and per-SNP sample sizes (equivalent to the eQTLGen-reported NrSamples; see results/M3_recalibrated_IV_table.csv).
- The transcriptome DE uses the deposited normalised counts (linear scale) with a 0.01 pseudo-count for
  log2 fold changes; see 02_gse289843_de.py.
- rs114796149 is absent from the FinnGen R5 release; the recalibrated PKD2 whole-blood analysis retains
  two instruments.
- The clinical-layer instrument files (`PKD_Clinical_22_IVs_Details.csv`, `outcome_dat_22IV_combined.csv`)
  store the harmonisation-stage snapshot of 22 SNPs (two of the 24 extracted SNPs were already excluded
  for allele mismatches with the FinnGen R5 reference). rs2541143 is excluded inside the analysis scripts
  as a palindromic variant with intermediate allele frequency, yielding the 21-IV analysis set reported
  in the manuscript (see the script output `dropped: rs2541143`).

## Independent validation (CPTAC-CCRCC)

`transcriptome/05_cptac_validation.py` reproduces the manuscript's validation in the CPTAC
clear-cell renal cell carcinoma proteogenomic cohort (110 tumours vs 84 normal tissues at the
protein level; 110 vs 75 at the RNA level; no sample overlap with TCGA-KIRC):

- Direction-of-effect concordance for the 43-gene panel: 30/43 significant at protein level
  (BH FDR < 0.05), 36/43 at RNA level; core findings (SETD2/PBRM1/VHL downregulation,
  HAVCR1 upregulation, NPHS1/NPHS2 renal-identity loss) concordant at both layers.
- STING1 and TBK1 show protein-level upregulation (FDR = 1.3e-29 / 1.6e-26), supporting
  cytosolic DNA sensing activation.
- Mutation-expression correspondence: tumours with inactivating BAP1 (n=17), SETD2 (n=13) or
  PBRM1 (n=47) mutations have significantly lower own-gene protein abundance than wild-type
  (panel-level BH FDR < 4.9e-4), consistent with loss-of-function.

Inputs are downloaded automatically by the script (LinkedOmics CPTAC-CCRCC exports and the
CPTAC pan-cancer clinical metadata file).

## License

MIT

## Outputs

Script outputs (e.g., S5_sensitivity.tsv) are generated on execution; precomputed result files are provided under results/.
