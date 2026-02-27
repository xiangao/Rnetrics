#' Logit estimation with robust or cluster-robust standard errors
#'
#' Computes maximum likelihood estimates for the logistic regression model
#' \code{Pr(D=1|X) = exp(X delta) / (1 + exp(X delta))} and returns the
#' score matrix needed for dyadic-robust variance calculations. Port of
#' Graham's Python \code{ipt.logit}.
#'
#' @param D       Binary numeric vector (length n): outcome variable.
#' @param X       n x K numeric matrix or data frame: regressors. A constant is
#'   prepended unless \code{nocons = TRUE}.
#' @param c_id    Optional vector (length n): cluster identifiers for
#'   cluster-robust standard errors.
#' @param s_wgt   Optional numeric vector (length n): sampling weights.
#' @param nocons  Logical. If \code{TRUE}, do NOT prepend a constant.
#' @param silent  Logical. If \code{TRUE}, suppress printed output.
#' @param full    Logical. If \code{TRUE} (default), print coefficient table.
#'
#' @return A list with components:
#'   \describe{
#'     \item{beta}{K x 1 matrix of ML estimates.}
#'     \item{vcov}{K x K variance-covariance matrix.}
#'     \item{hess_logl}{K x K negative Hessian of log-likelihood.}
#'     \item{score_i}{n x K matrix of score contributions.}
#'     \item{ehat}{n x 1 matrix of fitted probabilities Pr(D=1|X).}
#'     \item{converged}{Logical: did optimisation converge?}
#'   }
#' @importFrom stats optim
#' @export
logit_fit <- function(D, X, c_id = NULL, s_wgt = NULL, nocons = FALSE,
                      silent = FALSE, full = TRUE) {

  # --- Coerce to matrix -------------------------------------------------------
  if (is.data.frame(X)) {
    var_names <- names(X)
    X <- as.matrix(X)
  } else {
    X <- as.matrix(X)
    var_names <- colnames(X)
    if (is.null(var_names)) var_names <- paste0("X", seq_len(ncol(X)))
  }

  dep_var <- if (!is.null(names(D))) names(D)[1L] else "D"
  D <- as.numeric(D)
  n <- length(D)
  K <- ncol(X)

  # --- Sampling weights -------------------------------------------------------
  if (is.null(s_wgt)) {
    sw <- rep(1.0, n)
  } else {
    sw <- as.numeric(s_wgt) / mean(s_wgt)
  }

  # --- Add constant -----------------------------------------------------------
  if (!nocons) {
    X <- cbind(Constant = 1.0, X)
    var_names <- c("Constant", var_names)
    K <- K + 1L
  }

  # --- Log-likelihood and gradient functions ----------------------------------
  logit_nll <- function(delta) {
    Xd  <- as.numeric(X %*% delta)
    # numerically stable: log(1 + exp(x)) = x + log(1 + exp(-x)) for x>0
    lp1 <- ifelse(Xd >= 0, Xd + log1p(exp(-Xd)), log1p(exp(Xd)))
    -sum(sw * (D * Xd - lp1))
  }

  logit_grad <- function(delta) {
    Xd     <- as.numeric(X %*% delta)
    p_hat  <- 1 / (1 + exp(-Xd))
    -as.numeric(crossprod(X, sw * (D - p_hat)))
  }

  logit_hess_fn <- function(delta) {
    Xd    <- as.numeric(X %*% delta)
    p_hat <- 1 / (1 + exp(-Xd))
    w     <- sw * p_hat * (1 - p_hat)
    crossprod(X * w, X)
  }

  # --- Starting values --------------------------------------------------------
  delta_sv <- rep(0.0, K)
  if (!nocons) {
    # Not reached when nocons=TRUE
    p0 <- mean(D)
    p0 <- max(1e-6, min(1 - 1e-6, p0))
    delta_sv[1L] <- log(p0 / (1 - p0))
  }

  # --- Optimisation (L-BFGS-B via optim, then refine with Newton) -------------
  ctrl <- if (silent) {
    list(maxit = 1000L, trace = 0L)
  } else {
    list(maxit = 10000L, trace = 2L)
  }

  res <- optim(
    par     = delta_sv,
    fn      = logit_nll,
    gr      = logit_grad,
    method  = "BFGS",
    control = ctrl
  )

  delta_ml  <- res$par
  converged <- (res$convergence == 0L)

  # --- Scores and Hessian at MLE ----------------------------------------------
  Xd        <- as.numeric(X %*% delta_ml)
  p_hat     <- 1 / (1 + exp(-Xd))
  ehat      <- matrix(p_hat, ncol = 1L)
  score_i   <- X * (sw * (D - p_hat))     # n x K (note: NOT negated — same sign as gradient)
  hess_logl <- logit_hess_fn(delta_ml)    # positive definite (info matrix sign)

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
    cat("-                     LOGIT ESTIMATION RESULTS                        -\n")
    cat("-----------------------------------------------------------------------\n")
    cat(sprintf("Dependent variable:        %s\n", dep_var))
    cat(sprintf("Number of observations, n: %d\n\n", n))
    print_coef(delta_ml, vcov_beta, var_names = var_names, alpha = 0.05)
    if (is.null(c_id)) {
      cat("NOTE: Huber-robust standard errors reported.\n")
    } else {
      cat("NOTE: Cluster-Huber-robust standard errors reported.\n")
      cat(sprintf("      Number of clusters = %d\n", NC))
    }
    if (!is.null(s_wgt)) cat("NOTE: (Sampling) Weighted MLE estimates computed.\n")
  }

  list(
    beta      = matrix(delta_ml, ncol = 1L),
    vcov      = vcov_beta,
    hess_logl = hess_logl,
    score_i   = score_i,
    ehat      = ehat,
    converged = converged
  )
}
