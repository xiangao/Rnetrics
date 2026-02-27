#' Poisson (QMLE) estimation with robust or cluster-robust standard errors
#'
#' Computes quasi-maximum likelihood estimates for the Poisson regression model
#' \code{E[Y|X] = exp(X phi)} and returns the score matrix needed for
#' dyadic-robust variance calculations. Port of Graham's Python
#' \code{ipt.poisson}.
#'
#' @param Y       Numeric vector (length n): count outcome variable.
#' @param X       n x K numeric matrix or data frame: regressors. A constant is
#'   prepended unless \code{nocons = TRUE}.
#' @param c_id    Optional vector (length n): cluster identifiers for
#'   cluster-robust standard errors.
#' @param s_wgt   Optional numeric vector (length n): sampling weights.
#' @param nocons  Logical. If \code{TRUE}, do NOT prepend a constant.
#' @param silent  Logical. If \code{TRUE}, suppress printed output.
#' @param full    Logical. If \code{TRUE} (default), print coefficient table.
#' @param phi_sv  Optional numeric vector of length K: starting values.
#'
#' @return A list with components:
#'   \describe{
#'     \item{beta}{K x 1 matrix of QMLE estimates.}
#'     \item{vcov}{K x K variance-covariance matrix.}
#'     \item{hess_logl}{K x K negative Hessian of log-likelihood.}
#'     \item{score_i}{n x K matrix of score contributions.}
#'     \item{ehat}{n x 1 matrix of fitted conditional means exp(X phi).}
#'     \item{converged}{Logical: did optimisation converge?}
#'   }
#' @importFrom stats optim
#' @export
poisson_fit <- function(Y, X, c_id = NULL, s_wgt = NULL, nocons = FALSE,
                        silent = FALSE, full = TRUE, phi_sv = NULL) {

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

  # --- Negative log-likelihood and gradient -----------------------------------
  poisson_nll <- function(phi) {
    Xphi <- as.numeric(X %*% phi)
    mu   <- exp(Xphi)
    -sum(sw * (Y * Xphi - mu))
  }

  poisson_grad <- function(phi) {
    Xphi <- as.numeric(X %*% phi)
    mu   <- exp(Xphi)
    -as.numeric(crossprod(X, sw * (Y - mu)))
  }

  poisson_hess_fn <- function(phi) {
    Xphi <- as.numeric(X %*% phi)
    mu   <- exp(Xphi)
    crossprod(X * (sw * mu), X)
  }

  # --- Starting values --------------------------------------------------------
  if (is.null(phi_sv)) phi_sv <- rep(0.0, K)

  # --- Optimisation -----------------------------------------------------------
  ctrl <- if (silent) {
    list(maxit = 1000L, trace = 0L)
  } else {
    list(maxit = 10000L, trace = 2L)
  }

  res <- optim(
    par     = phi_sv,
    fn      = poisson_nll,
    gr      = poisson_grad,
    method  = "BFGS",
    control = ctrl
  )

  phi_ml    <- res$par
  converged <- (res$convergence == 0L)

  # --- Scores and Hessian at MLE ----------------------------------------------
  Xphi      <- as.numeric(X %*% phi_ml)
  mu_hat    <- exp(Xphi)
  ehat      <- matrix(mu_hat, ncol = 1L)
  score_i   <- X * (sw * (Y - mu_hat))    # n x K
  hess_logl <- poisson_hess_fn(phi_ml)    # K x K, positive definite

  # --- Variance-covariance matrix ---------------------------------------------
  if (is.null(c_id)) {
    fsc   <- n / (n - K)
    omega <- fsc * crossprod(score_i)
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
  }

  iH        <- solve(hess_logl)
  vcov_beta <- iH %*% omega %*% t(iH)

  # --- Output -----------------------------------------------------------------
  if (full) {
    cat("\n")
    cat("-----------------------------------------------------------------------\n")
    cat("-                    POISSON ESTIMATION RESULTS                       -\n")
    cat("-----------------------------------------------------------------------\n")
    cat(sprintf("Dependent variable:        %s\n", dep_var))
    cat(sprintf("Number of observations, n: %d\n\n", n))
    print_coef(phi_ml, vcov_beta, var_names = var_names, alpha = 0.05)
    if (is.null(c_id)) {
      cat("NOTE: Huber-robust standard errors reported.\n")
    } else {
      cat("NOTE: Cluster-Huber-robust standard errors reported.\n")
      cat(sprintf("      Number of clusters = %d\n", NC))
    }
    if (!is.null(s_wgt)) cat("NOTE: (Sampling) Weighted QMLE estimates computed.\n")
  }

  list(
    beta      = matrix(phi_ml, ncol = 1L),
    vcov      = vcov_beta,
    hess_logl = hess_logl,
    score_i   = score_i,
    ehat      = ehat,
    converged = converged
  )
}
