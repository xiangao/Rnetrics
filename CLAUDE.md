# Rnetrics — project notes for Claude

R port of Bryan Graham's Python `ipt` and `netrics` packages. Implements
dyadic regression — OLS, logit, or Poisson conditional-mean models on
pairwise (dyadic) data — with the bias-corrected dyadic-robust standard
errors (`DR_bc`) of Graham (forthcoming, *Handbook of Econometrics*), which
build on the dyadic-robust jackknife variance of Aronow, Samii & Assenova
(2015, *Political Analysis*).

## What's where

- `R/dyadic_regression.R` — main entry point `dyadic_regression()`. Fits the
  point estimates via one of the three `*_fit()` helpers below, then builds
  the variance matrix from two pieces: `Sigma2` (raw score covariance) and
  `Sigma1` (Hajek-projection covariance, computed via chunked scatter-add so
  the full n x K score matrix never needs to be materialized). `cov =`
  selects among `"ind"` (independence/HC), `"DR"` (uncorrected dyadic-robust),
  `"DR_bc"` (bias-corrected, default). Handles both undirected
  (n = N(N-1)/2) and directed (n = N(N-1)) dyad tables — the directed case
  symmetrizes scores across reverse dyads before building `Sigma2`.
- `R/ols.R`, `R/logit.R`, `R/poisson.R` — standalone point-estimation
  routines (`ols_fit`, `logit_fit`, `poisson_fit`), each returning
  `beta`, `vcov`, `hess_logl`, `score_i`, `ehat`. These are ports of
  Graham's Python `ipt.ols` / `ipt.logit` / `ipt.poisson` and are also
  exported for standalone use with HC or cluster-robust SEs (independent of
  the dyadic machinery). logit/poisson use `optim(method = "BFGS")`.
- `R/print_coef.R` — formatted coefficient table, port of `ipt.print_coef`.
- `tests/testthat/test-dyadic_regression.R` — simulated undirected dyadic
  DGPs; checks OLS beta against `lm()`, PSD-ness of the `DR_bc` vcov, and
  correct N/n bookkeping. (One test comment references "CLAUDE.md known
  results" — that refers to a template inherited from a sibling project;
  this repo has no such fixture, ignore that comment.)
- `tests/testthat/test-ols.R` — `ols_fit` against `lm()` and
  `sandwich::vcovHC`, plus cluster-robust and score-matrix-shape checks.
- `verification/` — one-off scripts (`run_netrics.py`, `run_Rnetrics.R`)
  and output logs from comparing this package against Graham's original
  Python `netrics`, not part of the test suite. Per `README.md`: max
  |coef diff| 2.3e-15 (simulated OLS) and 1.8e-03 (Poisson gravity, driven
  by optimizer tolerance differences between R's `optim(BFGS)` and Python's
  `fmin_bfgs`, not a structural bug); max |SE diff| 8.8e-05.

## Relationship to `netrics-fast`

Confirmed related, not a coincidence of naming: `~/projects/software/netrics-fast`
is a **Python** reimplementation of the *same* estimator — same
`dyadic_regression(Y, R, id_i, id_j, directed=, cov=)` signature, same
chunked O(nK) scatter-add Hajek projection, same `DR_bc` default, same
Graham (forthcoming) citation. Rnetrics is the R-native counterpart; treat
algorithmic changes here as candidates to port there too (and vice versa)
since they are meant to stay numerically interchangeable.

## Running tests

```r
devtools::load_all(quiet = TRUE)
testthat::test_dir("tests/testthat")
```

Fast (well under a minute). Last run: 16/16 passing, 0 failures/warnings.

## Gaps

- No pkgdown site yet — only `README.md` + roxygen `man/*.Rd`. If asked to
  add docs infrastructure, set one up (see `feedback_package_docs_default`
  convention used elsewhere in this portfolio: deploy pkgdown, link it from
  the README instead of GitHub blob links).
- No vignette; `README.md` usage example is the only worked walkthrough.
