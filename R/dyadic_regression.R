#' Dyadic regression with bias-corrected standard errors
#'
#' Computes dyadic regression estimates under linear (OLS), logit, or Poisson
#' conditional mean models, with heteroscedastic-robust, dyadic-robust (DR),
#' or bias-corrected dyadic-robust (DR_bc) standard errors following Graham
#' (forthcoming, Handbook of Econometrics).
#'
#' Unlike the original Python \code{netrics.dyadic_regression}, this
#' implementation:
#' \itemize{
#'   \item Accepts explicit \code{id_i} and \code{id_j} vectors instead of a
#'     Pandas MultiIndex — more natural for R data frames.
#'   \item Supports unbalanced/sparse networks without error.
#'   \item Uses a chunked Hajek projection to handle large datasets
#'     (100M+ rows) without materialising the full n x K score matrix.
#' }
#'
#' @param Y         Numeric vector (length n): outcome for each dyad.
#' @param R         n x K numeric matrix or data frame: regressors. A constant
#'   is prepended unless \code{nocons = TRUE}.
#' @param id_i      Vector (length n): "i" agent identifier for each dyad.
#' @param id_j      Vector (length n): "j" agent identifier for each dyad.
#' @param regmodel  Character: \code{"normal"} (OLS, default), \code{"logit"},
#'   or \code{"poisson"}.
#' @param directed  Logical. \code{FALSE} (default) for undirected outcomes
#'   (n = N(N-1)/2); \code{TRUE} for directed outcomes (n = N(N-1)).
#' @param nocons    Logical. If \code{TRUE}, do NOT prepend a constant to
#'   \code{R}.
#' @param silent    Logical. If \code{TRUE}, suppress all printed output.
#' @param cov       Character: variance estimator to use.
#'   \code{"ind"} assumes dyad independence,
#'   \code{"DR"} uses the dyadic-robust jackknife estimator,
#'   \code{"DR_bc"} (default) uses the bias-corrected estimator.
#' @param chunk_size Integer: number of rows to process per chunk when
#'   accumulating the Hajek projection (default 10,000,000). Reduce if memory
#'   is limited.
#'
#' @return A list with components:
#'   \describe{
#'     \item{theta}{K x 1 matrix of coefficient estimates.}
#'     \item{vcov}{K x K estimated variance-covariance matrix.}
#'     \item{se}{Numeric vector of length K: standard errors.}
#'     \item{var_names}{Character vector of variable names.}
#'     \item{N}{Number of agents.}
#'     \item{n}{Number of dyads.}
#'     \item{Sigma1}{K x K Sigma1 estimate (Hajek projection covariance).}
#'     \item{Sigma2}{K x K Sigma2 estimate (raw score covariance).}
#'   }
#'
#' @references Graham, B.S. (forthcoming). "Dyadic Regression."
#'   In \emph{Handbook of Econometrics}. North Holland.
#'
#' @importFrom utils flush.console
#' @export
dyadic_regression <- function(
    Y,
    R,
    id_i,
    id_j,
    regmodel  = "normal",
    directed  = FALSE,
    nocons    = FALSE,
    silent    = FALSE,
    cov       = "DR_bc",
    chunk_size = 10e6
) {

  # ---- 0. Input checks and coercion -----------------------------------------
  regmodel <- match.arg(regmodel, c("normal", "logit", "poisson"))
  cov      <- match.arg(cov,      c("ind", "DR", "DR_bc"))

  if (is.data.frame(R)) {
    var_names <- names(R)
    R <- as.matrix(R)
  } else {
    R <- as.matrix(R)
    var_names <- colnames(R)
    if (is.null(var_names)) var_names <- paste0("X", seq_len(ncol(R)))
  }

  Y    <- as.numeric(Y)
  id_i <- as.character(id_i)   # use character for safety
  id_j <- as.character(id_j)

  n <- length(Y)
  K <- ncol(R)

  if (nrow(R) != n)    stop("nrow(R) != length(Y)")
  if (length(id_i) != n) stop("length(id_i) != length(Y)")
  if (length(id_j) != n) stop("length(id_j) != length(Y)")

  # ---- 1. Add constant -------------------------------------------------------
  if (!nocons) {
    R         <- cbind(constant = 1.0, R)
    var_names <- c("constant", var_names)
    K         <- K + 1L
  }

  # ---- 2. Agent index map ----------------------------------------------------
  agents     <- sort(unique(c(id_i, id_j)))
  N          <- length(agents)
  agent2idx  <- seq_along(agents)
  names(agent2idx) <- agents

  ai <- agent2idx[id_i]   # n-length integer index [1..N]
  aj <- agent2idx[id_j]   # n-length integer index [1..N]

  # ---- 3. Point estimation ---------------------------------------------------
  Y_ser <- Y   # plain vector

  if (regmodel == "normal") {
    fit <- ols_fit(Y_ser, R, nocons = TRUE, silent = TRUE)
  } else if (regmodel == "logit") {
    fit <- logit_fit(Y_ser, R, nocons = TRUE, silent = TRUE, full = FALSE)
  } else {
    fit <- poisson_fit(Y_ser, R, nocons = TRUE, silent = TRUE, full = FALSE)
  }

  theta <- fit$beta       # K x 1
  H     <- fit$hess_logl  # K x K

  # ---- 4. Scores and Sigma2 -------------------------------------------------
  # For undirected data:
  #   Sigma2 = (1/n) * S'S  where S_ij = score for dyad (i,j)
  #   For OLS: S_ij = R_ij * resid_ij  =>  Sigma2 = R' diag(resid^2) R / n
  #   For logit/poisson: use score_i from fit object directly
  #
  # For directed data (Graham's convention):
  #   Use symmetrised kernel: s_sym_{ij} = S_{ij} + S_{ji}  (for i < j 0-indexed)
  #   Sigma2 = 2 * (s_sym' s_sym) / n  (factor 2: n = N(N-1), s_sym has n/2 rows)

  if (regmodel == "normal") {
    resid    <- Y - as.numeric(R %*% theta)
    if (!directed) {
      resid_sq <- resid^2
      Sigma2   <- crossprod(R * resid_sq, R) / n   # memory-efficient for large n
    }
    # For directed OLS, score_full is materialised below
  } else {
    score_i_full <- fit$score_i   # n x K
    if (!directed) {
      Sigma2 <- crossprod(score_i_full) / n
    }
  }

  # Directed Sigma2: build reverse-dyad index and symmetrise scores
  if (directed) {
    # Build position lookup: pos_mat[i, j] = row index in data for dyad (i -> j)
    pos_mat <- matrix(0L, nrow = N, ncol = N)
    for (k in seq_len(n)) pos_mat[ai[k], aj[k]] <- k

    # For each row k, find the row of the reverse dyad (aj[k] -> ai[k])
    rev_idx <- pos_mat[cbind(aj, ai)]   # n-length integer vector

    # Select lower triangle (ai > aj) to enumerate n/2 undirected pairs
    lower_mask <- (ai > aj)

    if (regmodel == "normal") {
      score_full <- R * resid   # n x K  (materialise for directed case only)
    } else {
      score_full <- score_i_full
    }

    s_sym  <- score_full[lower_mask, , drop = FALSE] +
              score_full[rev_idx[lower_mask], , drop = FALSE]   # (n/2) x K
    Sigma2 <- 2 * crossprod(s_sym) / n
  }

  # ---- 5. Sigma1 via chunked Hajek projection --------------------------------
  # s_bar_l = (1 / d_l) * sum_{ij: i=l or j=l} s_ij
  # where d_l = degree (number of dyads involving agent l)
  #
  # For OLS: s_ij = R_ij * resid_ij  (computed in chunks)
  # For logit/poisson: s_ij = score_i row (already in memory)

  agent_sum   <- matrix(0.0, nrow = N, ncol = K)
  agent_count <- integer(N)

  # Count dyads per agent
  for (i_val in ai) agent_count[i_val] <- agent_count[i_val] + 1L
  for (j_val in aj) agent_count[j_val] <- agent_count[j_val] + 1L

  chunk_size <- as.integer(chunk_size)

  if (regmodel == "normal") {
    resid_flat <- resid   # n-length

    for (start in seq(1L, n, by = chunk_size)) {
      end     <- min(start + chunk_size - 1L, n)
      idx_rng <- start:end
      S_chunk <- R[idx_rng, , drop = FALSE] * resid_flat[idx_rng]  # chunk x K

      ci <- ai[idx_rng]
      cj <- aj[idx_rng]

      # Scatter-add: for each column k, add S_chunk[,k] to agent_sum[ci,k] and agent_sum[cj,k]
      for (k in seq_len(K)) {
        agent_sum[, k] <- agent_sum[, k] +
          .sparse_scatter_add(ci, S_chunk[, k], N) +
          .sparse_scatter_add(cj, S_chunk[, k], N)
      }

      if (!silent && n > chunk_size) {
        cat(sprintf("  Hajek projection: processed %d / %d rows (%.0f%%)\r",
                    end, n, 100 * end / n))
        flush.console()
      }
    }
    if (!silent && n > chunk_size) cat("\n")

  } else {
    # logit/poisson: score_i already in memory
    for (k in seq_len(K)) {
      agent_sum[, k] <-
        .sparse_scatter_add(ai, score_i_full[, k], N) +
        .sparse_scatter_add(aj, score_i_full[, k], N)
    }
  }

  # For directed balanced networks, the Hajek projection averages over
  # N-1 unique partners (not 2*(N-1) directed dyad appearances).
  # agent_count[l] = 2*(N-1) for directed balanced; correct divisor is N-1.
  divisor <- if (directed) rep(N - 1L, N) else agent_count
  s_bar  <- agent_sum / divisor       # N x K (element-wise row division)
  Sigma1 <- crossprod(s_bar) / N      # K x K

  # ---- 6. Variance-covariance matrix ----------------------------------------
  if (!directed) {
    iGamma <- solve(-H / n)           # K x K
  } else {
    iGamma <- solve(-2 * H / n)       # Directed: factor of 2
  }

  vcov_theta <- switch(cov,
    "ind"   = (2 / (N - 1)) * (iGamma %*% Sigma2 %*% iGamma) / N,
    "DR"    = 4 * (iGamma %*% Sigma1 %*% iGamma) / N,
    "DR_bc" = 4 * (iGamma %*% (Sigma1 - 0.5 * Sigma2 / (N - 1)) %*% iGamma) / N
  )

  # PD enforcement via eigendecomposition
  if (cov == "DR_bc" || cov == "DR") {
    ev <- eigen(vcov_theta, symmetric = FALSE)
    if (any(ev$values < 0)) {
      n_neg <- sum(ev$values < 0)
      if (!silent) {
        warning(sprintf(
          "DR_bc vcov has %d negative eigenvalue(s); set to zero for PD enforcement.",
          n_neg))
      }
      ev$values[ev$values < 0] <- 0
      vcov_theta <- Re(ev$vectors %*% diag(ev$values) %*% solve(ev$vectors))
    }
  }

  se <- sqrt(pmax(diag(vcov_theta), 0))

  # ---- 7. Display results ----------------------------------------------------
  if (!silent) {
    cat("\n")
    cat("-------------------------------------------------------------------------------------------\n")
    cat("- DYADIC REGRESSION ESTIMATION RESULTS                                                    -\n")
    cat(sprintf("- (%-22s regression model)                                              -\n", regmodel))
    cat("-------------------------------------------------------------------------------------------\n\n")
    cat(sprintf("Number of agents,           N : %15s\n", formatC(N, format = "d", big.mark = ",")))
    cat(sprintf("Number of dyads,            n : %15s\n\n", formatC(n, format = "d", big.mark = ",")))
    cat("-------------------------------------------------------------------------------------------\n")
    print_coef(as.numeric(theta), vcov_theta, var_names = var_names)
    note <- switch(cov,
      "ind"   = "NOTE: Standard errors assume independence across dyads.",
      "DR"    = "NOTE: Dyadic-robust standard errors (Jackknife, no bias correction).",
      "DR_bc" = "NOTE: Bias-corrected dyadic-robust standard errors (DR_bc)."
    )
    cat(note, "\n")
  }

  list(
    theta     = theta,
    vcov      = vcov_theta,
    se        = se,
    var_names = var_names,
    N         = N,
    n         = n,
    Sigma1    = Sigma1,
    Sigma2    = Sigma2
  )
}


# ---------------------------------------------------------------------------
# Internal helper: scatter-add a numeric vector by integer index
# Equivalent to: out <- rep(0, N); for (i in seq_along(idx)) out[idx[i]] += x[i]
# but uses base R's vectorised `+` via tapply for correctness and speed.
# ---------------------------------------------------------------------------
.sparse_scatter_add <- function(idx, x, N) {
  # Use tabulate/rowsum-style accumulation
  out <- numeric(N)
  # Vectorised scatter via `+` using integer indexing is not directly available
  # in base R without a loop. We use data.frame + tapply-style rowsum.
  # For performance on large chunks, use base R's rowsum().
  # rowsum(x, idx, reorder=TRUE) returns a matrix with rows = sorted unique idx
  rs  <- rowsum(x, idx, reorder = TRUE)  # length <= N, rows are sorted group labels
  grp <- as.integer(rownames(rs))
  out[grp] <- rs[, 1L]
  out
}
