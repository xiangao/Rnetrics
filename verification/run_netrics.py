#!/usr/bin/env python
"""run_netrics.py

Run Graham's Python `netrics` package on two datasets:
  1. Simulated undirected dyadic data  (OLS, DR_bc)   — generated here, saved as CSV
  2. Log of Gravity (Santos Silva & Tenreyro 2006)    (Poisson, directed, DR_bc)

Saves coefficient / SE tables to CSV in this directory so they can be
loaded by run_Rnetrics.R for a side-by-side comparison.

Run from the project root:
  .venv/bin/python R/verification/run_netrics.py
"""
import sys, os

# ── locate ipt and netrics relative to project root ────────────────────────
ROOT   = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
OUTDIR = os.path.join(ROOT, 'R', 'verification')
sys.path.insert(0, os.path.join(ROOT, 'ipt'))
sys.path.insert(0, os.path.join(ROOT, 'netrics'))

import numpy as np
import pandas as pd
import netrics

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  SIMULATED DATA  ── N = 80 agents, undirected, OLS, DR_bc
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("PART 1: Simulated undirected dyadic data (OLS, DR_bc)")
print("=" * 70)

np.random.seed(123)
N = 80

# All undirected pairs  (i < j), 1-indexed
pairs = [(i, j) for i in range(1, N + 1) for j in range(i + 1, N + 1)]
id_i  = np.array([p[0] for p in pairs])
id_j  = np.array([p[1] for p in pairs])
n     = len(pairs)

a_cov = np.random.randn(N)
x_ij  = (a_cov[id_i - 1] - a_cov[id_j - 1]) ** 2
Y     = 1.0 + 0.5 * x_ij + np.random.normal(0.0, 0.5, n)

# Save shared dataset
sim_df = pd.DataFrame({'id_i': id_i, 'id_j': id_j, 'x': x_ij, 'Y': Y})
sim_csv = os.path.join(OUTDIR, 'sim_data.csv')
sim_df.to_csv(sim_csv, index=False)
print(f"  Simulated data saved to {sim_csv}")
print(f"  N = {N} agents,  n = {n} undirected dyads")

# Prepare MultiIndex for netrics
sim_mi = sim_df.set_index(['id_i', 'id_j']).sort_index()
Y_sim  = sim_mi['Y']
R_sim  = sim_mi[['x']]

[theta_s, vcov_s] = netrics.dyadic_regression(
    Y_sim, R_sim,
    regmodel='normal', directed=False, nocons=False, silent=True, cov='DR_bc'
)

se_s = np.sqrt(np.diag(vcov_s))
print("\n  Python netrics results (DR_bc):")
print(f"  {'Variable':<12}  {'Coef':>12}  {'SE':>12}")
print("  " + "-" * 40)
for name, c, s in zip(['constant', 'x'], theta_s.flatten(), se_s):
    print(f"  {name:<12}  {c:12.6f}  {s:12.6f}")

pd.DataFrame({
    'variable': ['constant', 'x'],
    'coef':      theta_s.flatten(),
    'se':        se_s
}).to_csv(os.path.join(OUTDIR, 'netrics_sim_results.csv'), index=False)


# ═══════════════════════════════════════════════════════════════════════════════
# 2.  LOG OF GRAVITY  ── N = 136 countries, directed, Poisson, DR_bc
#     Santos Silva & Tenreyro (2006, REStat 88(4): 641–658)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 2: Log of Gravity (Santos Silva & Tenreyro 2006)")
print("        Poisson, directed = True, DR_bc")
print("=" * 70)

grav_path = os.path.join(ROOT, 'short_courses', 'St_Gallen', '2024',
                         'Data', 'Log of Gravity.dta')
LogOfGravity = pd.read_stata(grav_path)
LogOfGravity.set_index(['s1_im', 's2_ex'], drop=False, inplace=True)
LogOfGravity.sort_index(level=['s1_im', 's2_ex'], inplace=True)

Y_grav = LogOfGravity['trade'].copy() / 1000   # scale as in Graham's notebook

COLS = ['lypex', 'lypim', 'lyex', 'lyim', 'ldist', 'border',
        'comlang', 'colony', 'landl_ex', 'landl_im',
        'lremot_ex', 'lremot_im', 'comfrt_wto', 'open_wto']
W_grav = LogOfGravity[COLS].copy()
W_grav['constant'] = 1.0   # constant last, as in Graham's notebook

N_grav = len(LogOfGravity.index.get_level_values('s1_im').unique())
n_grav = len(LogOfGravity)
print(f"  N = {N_grav} countries,  n = {n_grav} directed dyads")

[theta_g, vcov_g] = netrics.dyadic_regression(
    Y_grav, W_grav,
    regmodel='poisson', directed=True, nocons=True, silent=True, cov='DR_bc'
)

var_names_g = COLS + ['constant']
se_g = np.sqrt(np.diag(vcov_g))

print(f"\n  {'Variable':<15}  {'Coef':>12}  {'SE':>12}")
print("  " + "-" * 44)
for name, c, s in zip(var_names_g, theta_g.flatten(), se_g):
    print(f"  {name:<15}  {c:12.6f}  {s:12.6f}")

pd.DataFrame({
    'variable': var_names_g,
    'coef':      theta_g.flatten(),
    'se':        se_g
}).to_csv(os.path.join(OUTDIR, 'netrics_gravity_results.csv'), index=False)

print(f"\nResults written to {OUTDIR}/")
