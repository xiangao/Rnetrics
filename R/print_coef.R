#' Print coefficient table
#'
#' Prints a formatted table of coefficient estimates, standard errors, and
#' confidence intervals. Port of Graham's Python \code{ipt.print_coef}.
#'
#' @param beta   Numeric vector of length K: coefficient estimates.
#' @param vcov   K x K numeric matrix: estimated variance-covariance matrix.
#' @param var_names Character vector of length K: variable names. If \code{NULL},
#'   variables are labelled \code{X_0}, \code{X_1}, ...
#' @param alpha  Significance level for confidence interval (default 0.05 for 95\%).
#'
#' @return Invisibly returns \code{NULL}. All output is printed to the console.
#' @importFrom stats qnorm
#' @export
print_coef <- function(beta, vcov, var_names = NULL, alpha = 0.05) {
  K <- length(beta)

  if (is.null(var_names)) {
    var_names <- paste0("X_", seq(0, K - 1))
  }

  crit <- qnorm(1 - alpha / 2)
  se   <- sqrt(diag(as.matrix(vcov)))

  cat("\n")
  header <- sprintf(
    "%-25s %10s   (%10s)     (%s Confid. Interval )",
    "Independent variable", "Coef.", "Std. Err.", formatC(1 - alpha, digits = 2, format = "f")
  )
  cat(header, "\n")
  cat(strrep("-", 91), "\n")

  for (k in seq_len(K)) {
    lo <- beta[k] - crit * se[k]
    hi <- beta[k] + crit * se[k]
    cat(sprintf(
      "%-25s %10.6f   (%10.6f)     (%10.6f , %10.6f)\n",
      var_names[k], beta[k], se[k], lo, hi
    ))
  }

  cat("\n")
  cat(strrep("-", 91), "\n")
  invisible(NULL)
}
