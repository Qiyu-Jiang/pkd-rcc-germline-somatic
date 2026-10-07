# NOTE (updated 2026-09-15): the whole-blood eQTL layer was subsequently recalibrated
# from Z-scores using per-SNP allele frequencies (beta = Z / sqrt(2*p*(1-p)*N));
# see results/M3_recalibrated_IV_table.csv and transcriptome/04_eqtlgen_recalibration.py.
# This script documents the original preprocessing pipeline.
# ==============================================================================
# Supplementary Code — Mendelian randomization analysis
#
# Manuscript: Dissecting germline versus somatic mechanisms in the polycystic
#             kidney disease-renal cell carcinoma relationship: a multi-layer integrative study
#
# This single, self-contained script reproduces all Mendelian randomization
# results reported in the manuscript (Table 2, Table 3, Figs. 5-9) and the
# related supplementary tables (S8-S9).
#
# Cohort notes:
#   - TCGA-KIRC differential expression (533 tumours vs 72 adjacent normals)
#     was performed in Python and is NOT part of this script.
#   - The TCGA-KIRC 43-gene panel statistics in Table 1 correspond to the
#     533-tumour cohort (sample TCGA-DV-A4W0-05 excluded as a second primary
#     tumour; consistent with S1 and the genome-wide DE set).
#
# Requirements:
#   R >= 4.2
#   TwoSampleMR  >= 0.7.9  (github.com/MRCIEU/TwoSampleMR)
#   ieugwasr     >= 1.1.0
#   data.table, ggplot2
#   MRPRESSO 1.0   (github.com/MRCIEU/MRPresso) - optional, for Section 5
#
# Authentication:
#   Register at https://api.opengwas.io and set:
#     OPENGWAS_JWT <- "your full JWT here"
#   (Do NOT commit real tokens to the repository.)
#
# Data inputs (public; see manuscript Data Availability):
#   - Clinical PKD GWAS:            ebi-a-GCST90018904 (OpenGWAS)
#   - Renal cancer outcome:         finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC
#   - BMI positive control:         ukb-b-19953
#   - Whole-blood cis-eQTL:         eQTLGen consortium (2020)
#     https://www.eqtlgen.org/cis-eqtls.html
#   - GTEx v8 kidney cortex:
#     https://gtexportal.org/home/datasets
# ==============================================================================

## =============================================================================
## Section 0. Environment
## =============================================================================
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(ggplot2)

## Authentication (required for all OpenGWAS API calls)
OPENGWAS_JWT <- ""          # <- paste your full JWT here
if (!nzchar(OPENGWAS_JWT)) stop("Please set OPENGWAS_JWT (https://api.opengwas.io).")
Sys.setenv(OPENGWAS_JWT = OPENGWAS_JWT)

## Output directory (created if missing)
outdir <- "MR_output"
if (!dir.exists(outdir)) dir.create(outdir)

## =============================================================================
## Section 1. Clinical layer: IV extraction
##            Manuscript: Methods 2.4; Results 3.2.1
##            Exposure: clinical PKD diagnosis (ebi-a-GCST90018904; 424 cases)
##            Threshold P < 1e-5 (relaxed for the small exposure GWAS);
##            LD clumping r2 = 0.001, kb = 10,000 (EUR) -> 24 IVs
## =============================================================================
exposure_dat <- extract_instruments(
  outcomes = "ebi-a-GCST90018904",
  clump    = TRUE,
  p1       = 1e-5,
  r2       = 0.001,
  kb       = 10000
)
cat("IVs extracted:", nrow(exposure_dat), "(expected 24)\n")

## =============================================================================
## Section 2. Outcome extraction and harmonisation
##            3 IVs with no directly available outcome data are represented by
##            LD proxies (r2 > 0.8) as returned by the OpenGWAS API; 1 IV
##            (rs2541143) is dropped as palindromic with intermediate allele
##            frequency -> final 21 IVs (18 directly matched + 3 returned as LD proxies by the OpenGWAS API).
## =============================================================================
outcome_dat <- extract_outcome_data(
  snps     = exposure_dat$SNP,
  outcomes = "finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC"
)
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- subset(dat, mr_keep)
cat("IVs after harmonisation:", nrow(dat), "(expected 21)\n")

## =============================================================================
## Section 3. Primary MR analysis - Table 2
##            All estimators are computed within TwoSampleMR (a single
##            implementation); the multiplicative random-effects IVW is the
##            primary estimator owing to detected heterogeneity.
## =============================================================================
mr_res <- mr(dat)
res_table <- generate_odds_ratios(mr_res)
print(res_table)
fwrite(res_table, file.path(outdir, "Table2_clinical_MR_21IVs.csv"))

## =============================================================================
## Section 4. Heterogeneity and pleiotropy
## =============================================================================
het   <- mr_heterogeneity(dat);   print(het)
pleio <- mr_pleiotropy_test(dat); print(pleio)
fwrite(het,   file.path(outdir, "Table2_heterogeneity.csv"))
fwrite(pleio, file.path(outdir, "Table2_pleiotropy.csv"))

## =============================================================================
## Section 5. MR-PRESSO (18 directly matched IVs, proxies excluded)
## =============================================================================
presso_snps <- dat$SNP[dat$proxy.outcome == FALSE]
if (length(presso_snps) > 0 && requireNamespace("MRPRESSO", quietly = TRUE)) {
  presso_out <- MRPRESSO::mr_presso(
    BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome", SdExposure = "se.exposure",
    data = subset(dat, SNP %in% presso_snps),
    NbDistribution = 1000, SignifThreshold = 0.05)
  print(presso_out[["Global Test"]])
  fwrite(as.data.frame(presso_out[["Global Test"]]),
         file.path(outdir, "MRPresso_Global_18IVs.csv"))
}

## =============================================================================
## Section 6. Sensitivity analysis: 18 directly matched IVs
## =============================================================================
dat_sens <- subset(dat, !(SNP %in% dat$SNP[dat$proxy.outcome == TRUE]))
res19 <- generate_odds_ratios(mr(dat_sens))
print(res19)
fwrite(res19, file.path(outdir, "Table2_sensitivity_18IVs.csv"))

## =============================================================================
## Section 7. Leave-one-out analysis (clinical layer)
## =============================================================================
loo <- mr_leaveoneout(dat)
fwrite(loo, file.path(outdir, "LOO_clinical.csv"))
p1 <- mr_leaveoneout_plot(loo)[[1]]
ggsave(file.path(outdir, "Fig_LOO_clinical.pdf"), p1, width = 8, height = 6)
p2 <- mr_scatter_plot(mr_res, dat)[[1]]
ggsave(file.path(outdir, "Fig_scatter_clinical.pdf"), p2, width = 8, height = 6)

## =============================================================================
## Section 8. Power calculation (observed-precision approach)
##            Manuscript: Methods 2.8; Results 3.2.1
##            Replaces the earlier liability-scale R2 approach: power is
##            evaluated from the observed SE of the primary (multiplicative
##            random-effects) estimate.
## =============================================================================
i_re <- grepl("multiplicative random effects", res_table$method, ignore.case = TRUE)
se_primary <- res_table$se[i_re][1]
or80 <- exp((qnorm(0.975) + qnorm(0.80)) * se_primary)
cat("Observed SE (primary estimate):", round(se_primary, 4), "\n")
cat("80% power minimum detectable OR:", round(or80, 3), "\n")
# Manuscript wording: "the analysis had approximately 80% power to detect
# OR >= 1.16; effects of smaller magnitude cannot be excluded."

## Per-SNP instrument strength: F = (beta/SE)^2 (mean F = 21)
dat$F_per_snp <- (dat$beta.exposure / dat$se.exposure)^2
cat("Mean F-statistic:", round(mean(dat$F_per_snp), 2), "\n")
fwrite(dat[, c("SNP","beta.exposure","se.exposure","pval.exposure","F_per_snp")],
       file.path(outdir, "IV_characteristics_21IVs.csv"))

## =============================================================================
## Section 9. Whole-blood eQTL layer (PKD1 + PKD2 cis-eQTLs) - Table 3
##            Manuscript: Methods 2.5; Results 3.2.3
##   Input: PKD12_eQTL_subset.txt - locally filtered from the eQTLGen 2019-12-11
##          cis-eQTL file (FDR < 0.05, cis window, PKD1 ENSG00000008710 and
##          PKD2 ENSG00000118762). Columns (tab-separated, no header):
##          Pvalue, SNP, SNPChr, SNPPos, AssessedAllele, OtherAllele, Zscore,
##          Gene, GeneSymbol, GeneChr, GenePos, NrCohorts, NrSamples
##   Note (Limitations): SE = 1/sqrt(N) assumes MAF = 0.5; this underestimates
##   SE by 29-58% for MAF 0.5-0.1 and biases Wald ratios away from the null.
## =============================================================================
## =============================================================================
## REVISION NOTE (2026-09-15): the SE = 1/sqrt(N) approximation implemented in
## this section has been superseded. The whole-blood estimates reported in
## Table 3 were produced with per-SNP recalibration,
##   beta = Z / sqrt(2 * p * (1 - p) * N),  SE = 1 / sqrt(2 * p * (1 - p) * N),
## implemented in transcriptome/04_eqtlgen_recalibration.py
## (output: results/M3_recalibrated_IV_table.csv).
## The historical code below is disabled (wrapped in if (FALSE)); to run it,
## change `if (FALSE) {` to `if (TRUE) {`.
## =============================================================================
if (FALSE) {

eqtl <- fread("PKD12_eQTL_subset.txt", header = FALSE,
              col.names = c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele",
                            "OtherAllele","Zscore","Gene","GeneSymbol","GeneChr",
                            "GenePos","NrCohorts","NrSamples"))
eqtl_sig <- eqtl[Pvalue < 5e-8 & GeneSymbol %in% c("PKD1","PKD2")]
eqtl_sig[, SE_calc := 1/sqrt(as.numeric(NrSamples))]          # see Limitations
eqtl_sig[, Beta_calc := Zscore * SE_calc]

exp_dat <- format_data(as.data.frame(eqtl_sig), type = "exposure",
                       snp_col = "SNP", beta_col = "Beta_calc",
                       se_col = "SE_calc", pval_col = "Pvalue",
                       effect_allele_col = "AssessedAllele",
                       other_allele_col = "OtherAllele",
                       gene_col = "GeneSymbol", id_col = "GeneSymbol")
clumped <- ld_clump(data.frame(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure),
                    clump_r2 = 0.001, clump_kb = 10000, pop = "EUR")   # 5 IVs
exp_cl <- exp_dat[exp_dat$SNP %in% clumped$rsid, ]
out_blood <- extract_outcome_data(snps = exp_cl$SNP,
                                  outcomes = "finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC")
dat_blood <- harmonise_data(exp_cl, out_blood)
dat_blood <- subset(dat_blood, mr_keep)
cat("Blood-layer IVs after harmonisation:", nrow(dat_blood), "(expected 4)\n")

res_blood <- generate_odds_ratios(mr(dat_blood))
print(res_blood)
fwrite(res_blood, file.path(outdir, "Table3_blood_eQTL_MR.csv"))

## Per-gene analyses: PKD1 (single-IV Wald ratio) and PKD2 (3 IVs)
for (g in c("PKD1","PKD2")) {
  dg <- subset(dat_blood, grepl(g, exposure, ignore.case = TRUE))
  if (nrow(dg) > 0) {
    rg <- generate_odds_ratios(mr(dg))
    print(rg)
    fwrite(rg, file.path(outdir, paste0("Table3_", g, "_blood_eQTL.csv")))
  }
}
}
## =============================================================================
## Section 10. GTEx v8 kidney cortex: documented negative
##             Manuscript: Methods 2.4; Results 3.2.2
##   Input: GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.signif_variant_gene_pairs.txt.gz
##   (download from https://gtexportal.org/home/datasets)
## =============================================================================
gtex_file <- "Kidney_Cortex.v8.signif_variant_gene_pairs.txt.gz"
if (file.exists(gtex_file)) {
  gtex <- fread(gtex_file)
  target <- c("ENSG00000008710", "ENSG00000118762")   # PKD1, PKD2
  gtex_sig <- gtex[gene_id %in% target & pval_nominal < 5e-8]
  cat("GTEx v8 kidney-cortex significant cis-eQTLs for PKD1/PKD2:", nrow(gtex_sig), "\n")
  cat("(Expected 0; the kidney-cortex eQTL layer was therefore not feasible.)\n")
} else {
  cat("[Section 10 skipped] GTEx file not found in working directory.\n")
}

## =============================================================================
## Section 11. Positive control: BMI -> renal cancer (419 IVs)
##             Manuscript: Results 3.2.1; Table S9
## =============================================================================
bmi_exp <- extract_instruments(outcomes = "ukb-b-19953", clump = TRUE)
bmi_out <- extract_outcome_data(snps = bmi_exp$SNP,
                                outcomes = "finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC")
bmi_dat <- harmonise_data(bmi_exp, bmi_out)
bmi_res <- generate_odds_ratios(mr(bmi_dat))
print(generate_odds_ratios(mr(bmi_dat)))
fwrite(bmi_res, file.path(outdir, "Table_S9_positive_control_BMI.csv"))
# Expected: IVW OR = 1.42, 95% CI 1.06-1.90, P = 0.018
