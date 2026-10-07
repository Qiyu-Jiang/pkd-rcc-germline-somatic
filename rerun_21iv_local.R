# ==============================================================================
# rerun_21iv_local.R — clinical-layer MR on the final 21-IV set (manual allele harmonisation).
# Computes MR-Egger intercept, MR-PRESSO (reference) and leave-one-out statistics.
# ==============================================================================
suppressPackageStartupMessages({ library(TwoSampleMR); library(data.table) })

DD <- "./"   # working directory = folder containing the two CSV files
exp <- fread(paste0(DD, "PKD_Clinical_22_IVs_Details.csv"))
out <- fread(paste0(DD, "outcome_dat_22IV_combined.csv"))
m <- merge(exp, out, by = "SNP")
cat("merged:", nrow(m), "(expected: 22)\n")

m$beta.outcome.h <- ifelse(
  m$effect_allele.exposure == m$effect_allele.outcome, m$beta.outcome,
  ifelse(m$effect_allele.exposure == m$other_allele.outcome, -m$beta.outcome, NA_real_))

# rs2541143 excluded (palindromic, intermediate allele frequency).
d21 <- subset(m, SNP != "rs2541143" & !is.na(beta.outcome.h))
cat("dropped:", setdiff(m$SNP, d21$SNP), "\n")
cat("final IVs:", nrow(d21), "(expected: 21)\n")
stopifnot(nrow(d21) == 21)

dat <- data.frame(
  SNP = d21$SNP,
  beta.exposure = d21$beta.exposure, se.exposure = d21$se.exposure,
  effect_allele.exposure = d21$effect_allele.exposure,
  other_allele.exposure = d21$other_allele.exposure,
  beta.outcome = d21$beta.outcome.h, se.outcome = d21$se.outcome,
  effect_allele.outcome = d21$effect_allele.exposure,
  other_allele.outcome = d21$other_allele.exposure,
  pval.exposure = d21$pval.exposure,
  units.exposure = "log odds", units.outcome = "log odds",
  id.exposure = "PKD", id.outcome = "RCC",
  exposure = "clinical PKD", outcome = "renal cancer",
  stringsAsFactors = FALSE)
dat$mr_keep <- TRUE

## Sanity check: Q = 33.19 / df = 20 / P = 0.032 expected.
het <- mr_heterogeneity(dat)
print(het[, c("method", "Q", "Q_df", "Q_pval")])
write.csv(het, paste0(DD, "RERUN_het_21IV.csv"), row.names = FALSE)

## MR-Egger intercept.
pleio <- mr_pleiotropy_test(dat)
print(pleio)
write.csv(pleio, paste0(DD, "RERUN_pleiotropy_21IV.csv"), row.names = FALSE)

## MR-PRESSO (21 IVs, for reference; the sensitivity analysis uses the 18 directly matched IVs).
if (requireNamespace("MRPRESSO", quietly = TRUE)) {
  pr <- tryCatch(MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      data = dat, NbDistribution = 1000, SignifThreshold = 0.05),
    error = function(e) { cat("[PRESSO failed]", conditionMessage(e), "\n"); NULL })
  if (!is.null(pr)) {
    res <- pr[["MR Presso results"]]
    print(res[res$`MR method` %in% c("Global test", "Outlier-corrected"), ])
  }
}

## LOO
loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
if (!is.null(loo)) write.csv(loo, paste0(DD, "RERUN_LOO_21IV.csv"), row.names = FALSE)

## Primary analysis: five MR methods (consistency check against Table 2).
mr_res <- generate_odds_ratios(mr(dat))
print(mr_res[, c("method", "nsnp", "b", "se", "pval")])
cat("\n== done ==\n")
