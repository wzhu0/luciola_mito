################################################################################
#
# R Script: Convergence assessment for RevBayes MCMC output.
#
# Usage: Rscript 05_revbayes/scripts/assess_convergence.R <prefix>
#
# Example:
#   Rscript 05_revbayes/scripts/assess_convergence.R Luciola_mito_timetree
#   Rscript 05_revbayes/scripts/assess_convergence.R Luciola_mito_rooted_nonclock
#   Rscript 05_revbayes/scripts/assess_convergence.R Luciola_mito_unrooted_non-clock
#
# Authors: Wenjie Zhu
#
################################################################################

library(convenience)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide a filename prefix")
prefix <- args[1]
#prefix <- "Luciola_mito_timetree"
log_dir <- "05_revbayes/output"

log_files  <- paste0(log_dir, "/", prefix, "_run_", 1:4, ".log")

missing <- log_files[!file.exists(log_files)]
if (length(missing) > 0) stop("Missing files:\n", paste(missing, collapse = "\n"))

conv <- checkConvergence(
  list_files = log_files,
  format     = "revbayes",
  control    = makeControl(burnin = 0.2)
)

print(conv)

if (!is.null(conv$failed)) {
  cat("\nFailed parameters:\n")
  print(conv$failed)
} else {
  cat("\nAll parameters passed convergence assessment.\n")
}

q()
