#!/usr/bin/env Rscript
# run_testthat.R
# Run the full Rnetrics testthat suite and summarise results.
# Run from project root:  Rscript R/verification/run_testthat.R

library(devtools)
library(testthat)

OUTDIR <- "R/verification"

# ---- Run the test suite -------------------------------------------------------
cat("Running Rnetrics testthat suite ...\n\n")
results <- devtools::test("R/")

# ---- Write summary to file ---------------------------------------------------
out_file <- file.path(OUTDIR, "testthat_summary.txt")
con <- file(out_file, open = "wt")
sink(con, split = TRUE)

cat(strrep("=", 70), "\n")
cat("  RNETRICS — TESTTHAT SUITE RESULTS\n")
cat(strrep("=", 70), "\n\n")

cat("Test files:\n")
cat("  R/tests/testthat/test-ols.R             (4 tests)\n")
cat("  R/tests/testthat/test-dyadic_regression.R  (6 tests)\n\n")

cat("Test descriptions:\n\n")
cat("  test-ols.R:\n")
cat("    1. ols_fit matches lm() coefficients\n")
cat("    2. ols_fit HC vcov is close to sandwich::vcovHC (HC1)\n")
cat("    3. ols_fit cluster-robust SE works without error\n")
cat("    4. ols_fit returns score_i with correct dimensions\n\n")
cat("  test-dyadic_regression.R:\n")
cat("    1. dyadic_regression OLS beta matches lm() on simulated data\n")
cat("    2. DR_bc vcov is positive semi-definite\n")
cat("    3. SE ordering — DR_bc > ind (for typical correlated dyadic data)\n")
cat("    4. N and n are reported correctly\n")
cat("    5. cov = 'DR' (no bias correction) runs without error\n")
cat("    6. DR_bc SE > HC-robust SE for network data\n\n")

# Print testthat summary
print(results)

sink()
close(con)

cat(sprintf("\nSummary written to %s\n", out_file))
