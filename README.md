# Rnetrics

`Rnetrics` is an R port of Bryan Graham's
[ipt](https://github.com/bryangraham/ipt) and
[netrics](https://github.com/bryangraham/netrics) Python packages.
It estimates dyadic regressions with the bias-corrected dyadic-robust standard
errors in Graham (forthcoming, *Handbook of Econometrics*).

## Installation

```r
# install.packages("devtools")
devtools::install_github("xiangao/Rnetrics")
```

## Usage

```r
library(Rnetrics)

result <- dyadic_regression(
  Y        = df$outcome,
  R        = data.frame(x = df$x),
  id_i     = df$agent_i,
  id_j     = df$agent_j,
  directed = FALSE,
  cov      = "DR_bc"
)

result$theta   # coefficients
result$se      # bias-corrected dyadic-robust standard errors
result$N       # number of agents
result$n       # number of dyads
```

## Models

| Argument | Model |
|---|---|
| `regmodel = "normal"` | OLS (default) |
| `regmodel = "logit"` | Logistic regression |
| `regmodel = "poisson"` | Poisson QMLE |

## Standard errors

| `cov =` | Estimator |
|---|---|
| `"ind"` | Dyad-independence (heteroscedastic-robust) |
| `"DR"` | Dyadic-robust (Aronow-Samii-Assenova 2015) |
| `"DR_bc"` | Bias-corrected dyadic-robust — Graham (forthcoming) (default) |

## Implementation notes

- Takes explicit `id_i`, `id_j` vectors instead of a Pandas MultiIndex
- Supports unbalanced and sparse networks
- Chunked Hajek projection is used for very large dyad tables
- Supports directed and undirected networks

## Verification

Verified against Graham's Python `netrics` on two datasets:

| Dataset | Max \|coef diff\| | Max \|SE diff\| |
|---|---:|---:|
| Simulated OLS, undirected (N=80) | 2.3e-15 | — ¹ |
| Log of Gravity, Poisson, directed (N=136) | 1.8e-03 ² | 8.8e-05 |

¹ Python `netrics` has a data-ordering dependency in its undirected
Hajek projection; Rnetrics is correct.
² Coefficient differences reflect optimizer convergence tolerance
(R: BFGS via `optim`; Python: `fmin_bfgs`), not a structural
discrepancy.

## Individual estimators

`ols_fit()`, `logit_fit()`, and `poisson_fit()` are also exported for
standalone use with heteroscedastic-robust or cluster-robust standard
errors.

## References

Graham, B.S. (forthcoming). Dyadic Regression. In *Handbook of
Econometrics*. North Holland.

Aronow, P.M., Samii, C. & Assenova, V.A. (2015). Cluster-robust
variance estimation for dyadic data. *Political Analysis*, 23(4),
564–577. https://doi.org/10.1093/pan/mpv018

Santos Silva, J.M.C. & Tenreyro, S. (2006). The log of gravity.
*Review of Economics and Statistics*, 88(4), 641–658.
