#' OLS estimation with robust or cluster-robust standard errors
#'
#' Computes OLS coefficient estimates and returns the score matrix needed
#' for dyadic-robust variance calculations. Port of Graham's Python
#' \code{ipt.ols}.
#'
#' @param Y       Numeric vector (length n): dependent variable.
#' @param X       n x K numeric matrix or data frame: regressors. A constant is
#'   prepended unless \code{nocons = TRUE}.
#' @param c_id    Optional vector (length n): cluster identifiers for
#'   cluster-robust standard errors. If \code{NULL}, heteroscedastic-robust
#'   (HC) standard errors are reported.
#' @param s_wgt   Optional numeric vector (length n): sampling weights,
#'   normalised internally to have mean one.
#' @param nocons  Logical. If \code{TRUE}, do NOT prepend a constant column.
#'   Set to \code{TRUE} only when \code{X} already contains a constant.
#' @param silent  Logical. If \code{TRUE}, suppress printed output.
#'
#' @return A list with components:
#'   \describe{
#'     \item{beta}{K x 1 matrix of OLS estimates.}
#'     \item{vcov}{K x K variance-covariance matrix.}
#'     \item{hess_logl}{K x K Hessian (= X'X with weights).}
#'     \item{score_i}{n x K matrix of score contributions.}
#'     \item{ehat}{n x 1 matrix of fitted values.}
#'   }
#' @export
ols_fit <- function(Y, X, c_id = NULL, s_wgt = NULL, nocons = FALSE, silent = FALSE) {

  # --- Coerce to matrix -------------------------------------------------------
  if (is.data.frame(X)) {
    var_names <- names(X)
    X <- as.matrix(X)
  } else {
    X <- as.matrix(X)
    var_names <- colnames(X)
    if (is.null(var_names)) var_names <- paste0("X", seq_len(ncol(X)))
  }

  dep_var <- if (!is.null(names(Y))) names(Y)[1L] else "Y"
  Y <- as.numeric(Y)
  n <- length(Y)
  K <- ncol(X)

  # --- Sampling weights -------------------------------------------------------
  if (is.null(s_wgt)) {
    sw <- rep(1.0, n)
  } else {
    sw <- as.numeric(s_wgt) / mean(s_wgt)
  }

  # --- Add constant -----------------------------------------------------------
  if (!nocons) {
    X <- cbind(constant = 1.0, X)
    var_names <- c("constant", var_names)
    K <- K + 1L
  }

  # --- OLS point estimates ----------------------------------------------------
  swX  <- sw * X                          # n x K weighted X
  XX   <- crossprod(swX, X)              # K x K: (sw*X)' X
  XY   <- as.numeric(crossprod(swX, Y)) # K-vector: (sw*X)' Y
  beta <- solve(XX, XY)                  # K x 1
  ehat <- X %*% beta                     # n x 1

  # --- Scores -----------------------------------------------------------------
  resid   <- Y - as.numeric(ehat)
  score_i <- sw * X * resid              # n x K (broadcast)
  hess_logl <- XX

  # --- Variance-covariance matrix ---------------------------------------------
  if (is.null(c_id)) {
    fsc   <- n / (n - K)
    omega <- fsc * crossprod(score_i)
    iXX   <- solve(XX)
    vcov_beta <- iXX %*% omega %*% t(iXX)
  } else {
    c_list <- unique(c_id)
    NC     <- length(c_list)
    sum_score <- matrix(0.0, nrow = NC, ncol = K)
    for (ic in seq_len(NC)) {
      idx <- which(c_id == c_list[ic])
      sum_score[ic, ] <- colSums(score_i[idx, , drop = FALSE])
    }
    fsc   <- (n / (n - K)) * (NC / (NC - 1))
    omega <- fsc * crossprod(sum_score)
    iXX   <- solve(XX)
    vcov_beta <- iXX %*% omega %*% t(iXX)
  }

  # --- Output -----------------------------------------------------------------
  if (!silent) {
    cat("\n")
    cat("-----------------------------------------------------------------------\n")
    cat("-                     OLS ESTIMATION RESULTS                          -\n")
    cat("-----------------------------------------------------------------------\n")
    cat(sprintf("Dependent variable:        %s\n", dep_var))
    cat(sprintf("Number of observations, n: %d\n\n", n))
    print_coef(beta, vcov_beta, var_names = var_names, alpha = 0.05)
    if (is.null(c_id)) {
      cat("NOTE: Heteroscedastic-robust standard errors reported\n")
    } else {
      cat("NOTE: Cluster-robust standard errors reported\n")
      cat(sprintf("      Number of clusters = %d\n", NC))
    }
    if (!is.null(s_wgt)) {
      cat("NOTE: (Sampling) Weighted estimates computed.\n")
    }
  }

  list(
    beta      = matrix(beta, ncol = 1L),
    vcov      = vcov_beta,
    hess_logl = hess_logl,
    score_i   = score_i,
    ehat      = ehat
  )
}
