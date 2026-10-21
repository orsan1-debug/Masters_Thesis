# Verification of the stored summaries

Report written 2026-10-12 22:11 (R version 4.5.2 (2025-10-31 ucrt)). Each section lists the pairing checks run before storing, then the checks of verify.R. [x] passed; (info) never stops a build.

## correct_n (passed)

Built 2026-10-09 00:32; files: aipw_correct_n400k_dgp2.csv.gz, aipw_add_correct_n400k_dgp2.csv.gz.

- [x] add-on pairing (correct): oracle_hajek in aipw_add_correct_n400k_dgp2.csv.gz = aipw_correct_n400k_dgp2.csv.gz: 240000 of 240000 rows matched, max |diff| 0
- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 8 levels, 8 design rows
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 240 recomputed, 240 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.12, sbar_B 1.01 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 0.97, sbar_A 1.05, sbar_B 1.03 (B = 200)
- [x] (info) bootstrap/delta MCSE, aipw_glmnetcv vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 0.85, sbar_B 1.04 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs oracle_hajek, correct: bootstrap SD / delta MCSE: Rbar 1.05, sbar_A 1.01, sbar_B 1.10 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.00, sbar_B 0.99 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 80 recomputed, 80 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## misspec_n (passed)

Built 2026-10-09 00:32; files: aipw_misspec_n_dgp2.csv.gz, aipw_add_misspec_n_dgp2.csv.gz.

- [x] add-on pairing (misspec): oracle_hajek in aipw_add_misspec_n_dgp2.csv.gz = aipw_misspec_n_dgp2.csv.gz: 150000 of 150000 rows matched, max |diff| 1.6e-13
- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 5 levels, 5 design rows
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 180 recomputed, 180 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.02, sbar_A 0.97, sbar_B 1.06 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.96, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, aipw_glmnetcv vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.02, sbar_A 1.10, sbar_B 1.13 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs oracle_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 1.04, sbar_B 1.15 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs ols, misspec: bootstrap SD / delta MCSE: Rbar 0.99, sbar_A 1.05, sbar_B 1.03 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0_ortrue vs aipw_glmnetcv_ortrue, misspec: bootstrap SD / delta MCSE: Rbar 1.05, sbar_A 1.00, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, misspec: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 0.95, sbar_B 0.95 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 1.08, sbar_B 1.00 (B = 200)
- [ ] (info) bootstrap/delta MCSE, abw_balnet0 vs glmnetcv_ht, misspec: bootstrap SD / delta MCSE: Rbar 3.27, sbar_A 0.97, sbar_B 0.94 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 72 recomputed, 72 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## n_both (passed)

Built 2026-10-10 13:58; files: aipw_correct_n400k_dgp2.csv.gz, aipw_add_correct_n400k_dgp2.csv.gz, aipw_misspec_n_dgp2.csv.gz, aipw_add_misspec_n_dgp2.csv.gz.

- [x] add-on pairing (correct): oracle_hajek in aipw_add_correct_n400k_dgp2.csv.gz = aipw_correct_n400k_dgp2.csv.gz: 150000 of 150000 rows matched, max |diff| 0
- [x] add-on pairing (misspec): oracle_hajek in aipw_add_misspec_n_dgp2.csv.gz = aipw_misspec_n_dgp2.csv.gz: 150000 of 150000 rows matched, max |diff| 1.6e-13
- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 10 levels, 10 design rows
- [x] pooled specifications have the same levels and draws: 50000 and 50000 draws
- [x] shared draws across specifications: oracle_hajek: 150000 of 150000 rows matched, max |diff| 1.6e-13
- [x] shared draws across specifications: oracle_ht: 150000 of 150000 rows matched, max |diff| 1e-14
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 300 recomputed, 300 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 1.03, sbar_A 1.06, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 0.99, sbar_A 1.05, sbar_B 0.96 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, both: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.01, sbar_B 1.04 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.02, sbar_B 1.02 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 0.96, sbar_A 1.00, sbar_B 0.94 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 0.89, sbar_B 0.97 (B = 200)
- [x] (info) bootstrap/delta MCSE, aipw_glmnetcv vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 0.86, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, aipw_glmnetcv vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.08, sbar_A 0.94, sbar_B 0.98 (B = 200)
- [x] (info) bootstrap/delta MCSE, aipw_glmnetcv vs glmnetcv_hajek, both: bootstrap SD / delta MCSE: Rbar 0.84, sbar_A 0.93, sbar_B 0.95 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs oracle_hajek, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 0.99, sbar_B 1.15 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs oracle_hajek, misspec: bootstrap SD / delta MCSE: Rbar 0.99, sbar_A 0.89, sbar_B 0.94 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs oracle_hajek, both: bootstrap SD / delta MCSE: Rbar 0.96, sbar_A 0.89, sbar_B 1.09 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 1.02, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs ols, misspec: bootstrap SD / delta MCSE: Rbar 1.02, sbar_A 1.00, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs ols, both: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 1.03, sbar_B 1.01 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, correct: bootstrap SD / delta MCSE: Rbar 0.97, sbar_A 0.93, sbar_B 0.95 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, misspec: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.99, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, both: bootstrap SD / delta MCSE: Rbar 1.09, sbar_A 1.06, sbar_B 1.07 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 0.96, sbar_A 1.06, sbar_B 1.08 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.02, sbar_A 1.01, sbar_B 0.98 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs glmnetcv_hajek, both: bootstrap SD / delta MCSE: Rbar 1.03, sbar_A 1.02, sbar_B 1.10 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 0.88, sbar_A 0.95, sbar_B 0.95 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 0.99, sbar_A 0.90, sbar_B 1.04 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 0.99, sbar_A 0.94, sbar_B 1.01 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 192 recomputed, 192 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] both = correct + misspecified (cells): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## misspec_npath (passed)

Built 2026-10-10 10:42; files: aipw_full_misspec_nfine_dgp2.csv.gz, aipw_misspec_n_dgp2n1mill.csv.gz.

- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 30 levels, 30 design rows
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 630 recomputed, 630 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] (info) pairs: none in this section (per-cell measures only): 

## overlap (passed)

Built 2026-10-09 00:33; files: aipw_correct_overlap7_dgp2.csv.gz, aipw_correct_overlap8_dgp2.csv.gz, aipw_misspec_overlap7_dgp2.csv.gz, aipw_misspec_overlap8_dgp2.csv.gz.

- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 28 levels, 28 design rows
- [x] pooled specifications have the same levels and draws: 140000 and 140000 draws
- [x] shared draws across specifications: oracle_hajek: 420000 of 420000 rows matched, max |diff| 0
- [x] shared draws across specifications: oracle_ht: 420000 of 420000 rows matched, max |diff| 0
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 588 recomputed, 588 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 0.95, sbar_A 0.90, sbar_B 0.94 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 1.00, sbar_B 1.03 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.92, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 0.98, sbar_B 1.02 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 1.05, sbar_B 0.93 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, both: bootstrap SD / delta MCSE: Rbar 1.03, sbar_A 0.93, sbar_B 0.96 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, correct: bootstrap SD / delta MCSE: Rbar 1.09, sbar_A 1.05, sbar_B 1.07 (B = 200)
- [ ] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, misspec: bootstrap SD / delta MCSE: Rbar 4.45, sbar_A 1.00, sbar_B 1.02 (B = 200)
- [ ] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, both: bootstrap SD / delta MCSE: Rbar 3.23, sbar_A 1.12, sbar_B 0.96 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, correct: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 1.06, sbar_B 1.08 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, misspec: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.00, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, both: bootstrap SD / delta MCSE: Rbar 0.94, sbar_A 0.99, sbar_B 0.87 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 1.08, sbar_A 1.03, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.99, sbar_B 1.05 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.93, sbar_B 1.00 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 1200 recomputed, 1200 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] both = correct + misspecified (cells): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## nc (passed)

Built 2026-10-09 00:35; files: aipw_correct_n_by_overlap_dgp2.csv.gz, aipw_misspec_n_by_overlap_dgp2.csv.gz.

- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 1000: min 1000, max 1000
- [x] one design row per level: 312 levels, 312 design rows
- [x] pooled specifications have the same levels and draws: 156000 and 156000 draws
- [x] shared draws across specifications: oracle_hajek: 468000 of 468000 rows matched, max |diff| 0
- [x] shared draws across specifications: oracle_ht: 468000 of 468000 rows matched, max |diff| 0
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 5616 recomputed, 5616 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] (info) pairs: none in this section (per-cell measures only): 

## p_axis (passed)

Built 2026-10-09 10:21; files: aipw_correct_p_dgp2.csv.gz, aipw_misspec_p_dgp2.csv.gz.

- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 24 levels, 24 design rows
- [x] pooled specifications have the same levels and draws: 120000 and 120000 draws
- [x] shared draws across specifications: oracle_hajek: 360000 of 360000 rows matched, max |diff| 0
- [x] shared draws across specifications: oracle_ht: 360000 of 360000 rows matched, max |diff| 0
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 576 recomputed, 576 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] path ends: every cell found in both routes: 24 recomputed, 24 stored
- [x] path ends: means and MCSEs agree: largest |difference| between the routes
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 0.93, sbar_A 0.94, sbar_B 0.96 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.98, sbar_B 1.05 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 0.94, sbar_A 1.07, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.95, sbar_B 0.98 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 0.88, sbar_A 1.02, sbar_B 0.88 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, both: bootstrap SD / delta MCSE: Rbar 0.91, sbar_A 0.99, sbar_B 0.94 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, correct: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.97, sbar_B 0.96 (B = 200)
- [ ] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, misspec: bootstrap SD / delta MCSE: Rbar 1.98, sbar_A 1.02, sbar_B 1.92 (B = 200)
- [ ] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, both: bootstrap SD / delta MCSE: Rbar 1.68, sbar_A 1.03, sbar_B 1.81 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, correct: bootstrap SD / delta MCSE: Rbar 1.08, sbar_A 1.12, sbar_B 1.09 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, misspec: bootstrap SD / delta MCSE: Rbar 1.02, sbar_A 1.03, sbar_B 1.06 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, both: bootstrap SD / delta MCSE: Rbar 1.07, sbar_A 0.98, sbar_B 0.98 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 0.94, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 0.97, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 1.03, sbar_A 0.95, sbar_B 1.04 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 360 recomputed, 360 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] both = correct + misspecified (cells): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## sparsity (passed)

Built 2026-10-09 10:12; files: aipw_correct_sparsity_hd_dgp2.csv.gz, aipw_correct_sparsity_hd_equi_dgp3.csv.gz, aipw_misspec_tanobs_sparsity_hd_dgp3.csv.gz, aipw_misspec_tanobs_sparsity_hd_equi_dgp3.csv.gz.

- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 1000: min 1000, max 1000
- [x] one design row per level: 204 levels, 204 design rows
- [x] pooled specifications have the same levels and draws: 102000 and 102000 draws
- [x] shared draws across specifications: oracle_hajek: 204000 of 204000 rows matched, max |diff| 0
- [x] shared draws across specifications: oracle_ht: 204000 of 204000 rows matched, max |diff| 0
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 3264 recomputed, 3264 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] path ends: every cell found in both routes: 204 recomputed, 204 stored
- [x] path ends: means and MCSEs agree: largest |difference| between the routes
- [x] truncation: every cell found in both routes: 204 recomputed, 204 stored
- [x] truncation: shares and MCSEs agree: largest |difference| between the routes
- [x] cost: every stored row recomputed: 384 recomputed, 384 stored
- [x] cost: ratios, MCSEs, shares and cell counts agree: largest |difference| between the routes
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.08, sbar_B 1.01 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 0.84, sbar_A 1.00, sbar_B 0.95 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 1.06, sbar_A 1.03, sbar_B 1.05 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 0.96, sbar_A 0.99, sbar_B 0.98 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 0.87, sbar_A 1.00, sbar_B 0.96 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, both: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 1.02, sbar_B 1.02 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, correct: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 0.98, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, misspec: bootstrap SD / delta MCSE: Rbar 1.10, sbar_A 1.00, sbar_B 1.02 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, both: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.96, sbar_B 1.09 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, correct: bootstrap SD / delta MCSE: Rbar 1.10, sbar_A 1.08, sbar_B 1.10 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, misspec: bootstrap SD / delta MCSE: Rbar 1.07, sbar_A 1.17, sbar_B 1.12 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, both: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 0.98, sbar_B 0.98 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 0.94, sbar_B 0.96 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 0.97, sbar_B 0.97 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 0.95, sbar_A 1.11, sbar_B 1.12 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 4680 recomputed, 4680 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] both = correct + misspecified (cells): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## signal (passed)

Built 2026-10-09 02:16; files: aipw_correct_strength_dgp2.csv.gz, aipw_misspec_strength_dgp2.csv.gz.

- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 16 levels, 16 design rows
- [x] pooled specifications have the same levels and draws: 80000 and 80000 draws
- [x] shared draws across specifications: oracle_hajek: 240000 of 240000 rows matched, max |diff| 0
- [x] shared draws across specifications: oracle_ht: 240000 of 240000 rows matched, max |diff| 0
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 384 recomputed, 384 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 1.05, sbar_B 1.04 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 0.90, sbar_A 0.97, sbar_B 0.91 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 1.03, sbar_A 1.02, sbar_B 0.98 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.01, sbar_B 1.03 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.94, sbar_B 1.01 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, both: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 1.05, sbar_B 1.05 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, correct: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 1.01, sbar_B 0.96 (B = 200)
- [ ] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, misspec: bootstrap SD / delta MCSE: Rbar 1.36, sbar_A 0.98, sbar_B 1.57 (B = 200)
- [ ] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, both: bootstrap SD / delta MCSE: Rbar 1.37, sbar_A 0.96, sbar_B 1.63 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, correct: bootstrap SD / delta MCSE: Rbar 1.03, sbar_A 1.06, sbar_B 1.05 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, misspec: bootstrap SD / delta MCSE: Rbar 0.96, sbar_A 1.07, sbar_B 1.07 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, both: bootstrap SD / delta MCSE: Rbar 0.94, sbar_A 0.99, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 0.97, sbar_A 0.88, sbar_B 0.92 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 0.91, sbar_B 1.03 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 0.96, sbar_B 0.91 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 1320 recomputed, 1320 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] both = correct + misspecified (cells): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## treat_share (passed)

Built 2026-10-09 00:38; files: aipw_correct_treat_prop_dgp2.csv.gz, aipw_misspec_treat_prop_dgp2.csv.gz.

- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 16 levels, 16 design rows
- [x] pooled specifications have the same levels and draws: 80000 and 80000 draws
- [x] shared draws across specifications: oracle_hajek: 240000 of 240000 rows matched, max |diff| 0
- [x] shared draws across specifications: oracle_ht: 240000 of 240000 rows matched, max |diff| 0
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 384 recomputed, 384 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 1.02, sbar_A 1.02, sbar_B 0.97 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 1.10, sbar_A 1.00, sbar_B 1.09 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 1.07, sbar_A 1.06, sbar_B 1.03 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, correct: bootstrap SD / delta MCSE: Rbar 0.95, sbar_A 1.01, sbar_B 1.06 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, misspec: bootstrap SD / delta MCSE: Rbar 1.05, sbar_A 0.95, sbar_B 0.96 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_hajek, both: bootstrap SD / delta MCSE: Rbar 0.97, sbar_A 0.98, sbar_B 1.03 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, correct: bootstrap SD / delta MCSE: Rbar 0.92, sbar_A 1.00, sbar_B 0.93 (B = 200)
- [ ] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, misspec: bootstrap SD / delta MCSE: Rbar 2.58, sbar_A 1.02, sbar_B 2.33 (B = 200)
- [ ] (info) bootstrap/delta MCSE, balnet0 vs glmnetcv_ht, both: bootstrap SD / delta MCSE: Rbar 2.12, sbar_A 0.99, sbar_B 2.40 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, correct: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 1.11, sbar_B 1.08 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, misspec: bootstrap SD / delta MCSE: Rbar 0.98, sbar_A 1.00, sbar_B 1.03 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs balnet0, both: bootstrap SD / delta MCSE: Rbar 0.97, sbar_A 0.99, sbar_B 1.02 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, correct: bootstrap SD / delta MCSE: Rbar 1.10, sbar_A 1.00, sbar_B 1.05 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, misspec: bootstrap SD / delta MCSE: Rbar 1.02, sbar_A 1.07, sbar_B 1.04 (B = 200)
- [x] (info) bootstrap/delta MCSE, abw_balnet0 vs aipw_glmnetcv, both: bootstrap SD / delta MCSE: Rbar 0.94, sbar_A 1.04, sbar_B 1.01 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 120 recomputed, 120 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] both = correct + misspecified (cells): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## tol_n (passed)

Built 2026-10-09 00:39; files: aipw_correct_n400k_dgp2.csv.gz, aipw_add_correct_n400k_dgp2.csv.gz, aipw_tol_correct_n_dgp2.csv.gz.

- [x] add-on pairing (correct): oracle_hajek in aipw_add_correct_n400k_dgp2.csv.gz = aipw_correct_n400k_dgp2.csv.gz: 150000 of 150000 rows matched, max |diff| 0
- [x] add-on pairing (correct): oracle_hajek in aipw_tol_correct_n_dgp2.csv.gz = aipw_correct_n400k_dgp2.csv.gz: 150000 of 150000 rows matched, max |diff| 0
- [x] each arm read from one file only: no duplicated (spec, level, sim, outcome, estimator)
- [x] replications per arm and cell = 10000: min 10000, max 10000
- [x] one design row per level: 5 levels, 5 design rows
- [x] no infinite estimates: 0 infinite
- [x] perf: every arm-cell found in both routes: 180 recomputed, 180 stored
- [x] perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes
- [x] perf: missing counts agree: 
- [x] cells: each cell once per pair: 
- [x] (info) bootstrap/delta MCSE, balnet0_tol2 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.05, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol2.5 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 1.02, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol3 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.02, sbar_A 0.97, sbar_B 0.96 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol3.5 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.01, sbar_A 0.99, sbar_B 1.06 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.04, sbar_A 1.02, sbar_B 0.99 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol4.5 vs ols, correct: bootstrap SD / delta MCSE: Rbar 0.95, sbar_A 0.94, sbar_B 0.98 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol5 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.09, sbar_A 1.06, sbar_B 1.10 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol5.5 vs ols, correct: bootstrap SD / delta MCSE: Rbar 0.99, sbar_A 0.97, sbar_B 1.00 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol6 vs ols, correct: bootstrap SD / delta MCSE: Rbar 0.95, sbar_A 1.03, sbar_B 1.04 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol6.5 vs ols, correct: bootstrap SD / delta MCSE: Rbar 0.93, sbar_A 1.01, sbar_B 0.94 (B = 200)
- [x] (info) bootstrap/delta MCSE, balnet0_tol7 vs ols, correct: bootstrap SD / delta MCSE: Rbar 1.00, sbar_A 1.00, sbar_B 0.99 (B = 200)
- [x] cells: R_c, shares, RMSEs, MCSEs and kept flags agree: largest |difference| between the routes
- [x] averages: every stored row recomputed: 88 recomputed, 88 stored
- [x] averages: Rbar, sbar, MCSEs, cell counts agree: largest |difference| between the routes
- [x] pooled = sum over outcomes (cells and dropped): 
- [x] by-sets: every factor splits the same cells once: sum over a factor's levels of (cells + dropped) is the same for every factor of a set
- [x] common Theta: every pair averaged over the same number of cells: 
- [x] Rbar = 100 (1 - exp(mean log R_c)): 
- [x] (info) sbar of one arm across pairs (common Theta): differences come only from pairwise deletion of missing estimates

## app_misspec_n (passed)

Built 2026-10-09 00:40; files: aipw_misspec_n_dgp2.csv.gz, aipw_add_misspec_n_dgp2.csv.gz.

- [x] full (misspec): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (misspec): replications per arm and cell = 10000: min 10000, max 10000
- [x] diagnostics (misspec): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (misspec): replications per cell = 10000: min 10000, max 10000
- [x] approximation bias (misspec): abw_balnet0 and ols on the same draws, each cell once: 15 of 15 cells kept
- [x] full (misspec): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 135 arm-cells; largest |difference| between the routes
- [x] (info) full (misspec): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.01, middle 90% 0.90 to 1.13
- [x] diagnostics (misspec): medians, quartiles, SD, range agree: 160 rows; largest |difference| between the routes
- [x] diagnostics (misspec): path truncation counts agree: 
- [x] approximation bias (misspec): per-cell biases and MCSEs agree: 
- [x] approximation bias (misspec): averages over cells and MCSEs agree: 

## app_p_axis (passed)

Built 2026-10-09 00:41; files: aipw_correct_p_dgp2.csv.gz, aipw_misspec_p_dgp2.csv.gz.

- [x] full (correct): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (correct): replications per arm and cell = 10000: min 10000, max 10000
- [x] full (misspec): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (misspec): replications per arm and cell = 10000: min 10000, max 10000
- [x] diagnostics (correct): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (correct): replications per cell = 10000: min 10000, max 10000
- [x] diagnostics (misspec): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (misspec): replications per cell = 10000: min 10000, max 10000
- [x] full (correct): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 216 arm-cells; largest |difference| between the routes
- [x] (info) full (correct): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.01, middle 90% 0.91 to 1.17
- [x] full (misspec): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 288 arm-cells; largest |difference| between the routes
- [x] (info) full (misspec): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.01, middle 90% 0.89 to 1.15
- [x] diagnostics (correct): medians, quartiles, SD, range agree: 384 rows; largest |difference| between the routes
- [x] diagnostics (correct): path truncation counts agree: 
- [x] diagnostics (misspec): medians, quartiles, SD, range agree: 384 rows; largest |difference| between the routes
- [x] diagnostics (misspec): path truncation counts agree: 

## app_sparsity (passed)

Built 2026-10-09 00:43; files: aipw_correct_sparsity_hd_dgp2.csv.gz.

- [x] full (correct): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (correct): replications per arm and cell = 1000: min 1000, max 1000
- [x] diagnostics (correct): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (correct): replications per cell = 1000: min 1000, max 1000
- [x] full (correct): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 816 arm-cells; largest |difference| between the routes
- [x] (info) full (correct): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.00, middle 90% 0.89 to 1.13
- [x] diagnostics (correct): medians, quartiles, SD, range agree: 2040 rows; largest |difference| between the routes
- [x] diagnostics (correct): path truncation counts agree: 

## app_sparsity_equi (passed)

Built 2026-10-09 11:56; files: aipw_correct_sparsity_hd_equi_dgp3.csv.gz.

- [x] full (correct): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (correct): replications per arm and cell = 1000: min 1000, max 1000
- [x] diagnostics (correct): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (correct): replications per cell = 1000: min 1000, max 1000
- [x] full (correct): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 408 arm-cells; largest |difference| between the routes
- [x] (info) full (correct): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.00, middle 90% 0.89 to 1.16
- [x] diagnostics (correct): medians, quartiles, SD, range agree: 1020 rows; largest |difference| between the routes
- [x] diagnostics (correct): path truncation counts agree: 

## app_sparsity_tan (passed)

Built 2026-10-09 11:58; files: aipw_misspec_tanobs_sparsity_hd_dgp3.csv.gz, aipw_misspec_tanobs_sparsity_hd_equi_dgp3.csv.gz.

- [x] full (misspec): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (misspec): replications per arm and cell = 1000: min 1000, max 1000
- [x] diagnostics (misspec): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (misspec): replications per cell = 1000: min 1000, max 1000
- [x] full (misspec): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 1632 arm-cells; largest |difference| between the routes
- [x] (info) full (misspec): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.00, middle 90% 0.89 to 1.14
- [x] diagnostics (misspec): medians, quartiles, SD, range agree: 3060 rows; largest |difference| between the routes
- [x] diagnostics (misspec): path truncation counts agree: 

## app_nc (passed)

Built 2026-10-09 00:43; files: aipw_correct_n_by_overlap_dgp2.csv.gz, aipw_misspec_n_by_overlap_dgp2.csv.gz.

- [x] full (correct): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (correct): replications per arm and cell = 1000: min 1000, max 1000
- [x] full (misspec): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (misspec): replications per arm and cell = 1000: min 1000, max 1000
- [x] full (correct): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 2808 arm-cells; largest |difference| between the routes
- [x] full (misspec): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 3744 arm-cells; largest |difference| between the routes

## app_ipw_800k (passed)

Built 2026-10-09 06:13; files: correctspecDGP2.csv.gz, correctspecDGP400k.csv.gz, correctspecDGP800k.csv.gz.

- [x] full (correct): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (correct): replications per arm and cell = 1000: min 1000, max 1000
- [x] each draw read once: no duplicated (level, sim, arm, outcome)
- [x] replications per arm and cell = 1000: 36 rows; min 1000, max 1000 replications, 0 missing
- [x] full (correct): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 36 arm-cells; largest |difference| between the routes
- [x] (info) full (correct): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.01, middle 90% 0.89 to 1.16
- [x] simple: bias, EmpSE, RMSE and MCSEs agree: 36 arm-cells; largest |difference| between the routes

## app_correct_n (passed)

Built 2026-10-09 00:44; files: aipw_correct_n400k_dgp2.csv.gz, aipw_add_correct_n400k_dgp2.csv.gz.

- [x] full (correct): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (correct): replications per arm and cell = 10000: min 10000, max 10000
- [x] diagnostics (correct): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (correct): replications per cell = 10000: min 10000, max 10000
- [x] approximation bias (correct): abw_balnet0 and ols on the same draws, each cell once: 24 of 24 cells kept
- [x] full (correct): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 240 arm-cells; largest |difference| between the routes
- [x] (info) full (correct): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 0.99, middle 90% 0.88 to 1.13
- [x] diagnostics (correct): medians, quartiles, SD, range agree: 256 rows; largest |difference| between the routes
- [x] diagnostics (correct): path truncation counts agree: 
- [x] approximation bias (correct): per-cell biases and MCSEs agree: 
- [x] approximation bias (correct): averages over cells and MCSEs agree: 

## app_overlap (passed)

Built 2026-10-09 00:45; files: aipw_correct_overlap7_dgp2.csv.gz, aipw_diag_correct_overlap7_dgp2.csv.gz, aipw_misspec_overlap7_dgp2.csv.gz, aipw_full_misspec_overlap7_dgp2.csv.gz.

- [x] add-on pairing (correct): oracle_hajek in aipw_diag_correct_overlap7_dgp2.csv.gz = aipw_correct_overlap7_dgp2.csv.gz: 210000 of 210000 rows matched, max |diff| 0
- [x] approximation bias (correct): abw_balnet0 and ols on the same draws, each cell once: 21 of 21 cells kept
- [x] add-on pairing (misspec): oracle_hajek in aipw_full_misspec_overlap7_dgp2.csv.gz = aipw_misspec_overlap7_dgp2.csv.gz: 210000 of 210000 rows matched, max |diff| 0
- [x] approximation bias (misspec): abw_balnet0 and ols on the same draws, each cell once: 21 of 21 cells kept
- [x] approximation bias (correct): per-cell biases and MCSEs agree: 
- [x] approximation bias (correct): averages over cells and MCSEs agree: 
- [x] approximation bias (misspec): per-cell biases and MCSEs agree: 
- [x] approximation bias (misspec): averages over cells and MCSEs agree: 

## app_overlap_full (passed)

Built 2026-10-09 01:39; files: aipw_correct_overlap7_dgp2.csv.gz, aipw_correct_overlap8_dgp2.csv.gz, aipw_misspec_overlap7_dgp2.csv.gz, aipw_misspec_overlap8_dgp2.csv.gz.

- [x] full (correct): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (correct): replications per arm and cell = 10000: min 10000, max 10000
- [x] full (misspec): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (misspec): replications per arm and cell = 10000: min 10000, max 10000
- [x] diagnostics (correct): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (correct): replications per cell = 10000: min 10000, max 10000
- [x] diagnostics (misspec): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (misspec): replications per cell = 10000: min 10000, max 10000
- [x] full (correct): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 252 arm-cells; largest |difference| between the routes
- [x] (info) full (correct): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.00, middle 90% 0.89 to 1.15
- [x] full (misspec): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 252 arm-cells; largest |difference| between the routes
- [x] (info) full (misspec): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.00, middle 90% 0.89 to 1.14
- [x] diagnostics (correct): medians, quartiles, SD, range agree: 448 rows; largest |difference| between the routes
- [x] diagnostics (correct): path truncation counts agree: 
- [x] diagnostics (misspec): medians, quartiles, SD, range agree: 448 rows; largest |difference| between the routes
- [x] diagnostics (misspec): path truncation counts agree: 

## app_signal (passed)

Built 2026-10-09 02:18; files: aipw_correct_strength_dgp2.csv.gz, aipw_misspec_strength_dgp2.csv.gz.

- [x] full (correct): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (correct): replications per arm and cell = 10000: min 10000, max 10000
- [x] full (misspec): each arm read from one file: no duplicated (level, sim, arm)
- [x] full (misspec): replications per arm and cell = 10000: min 10000, max 10000
- [x] diagnostics (correct): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (correct): replications per cell = 10000: min 10000, max 10000
- [x] diagnostics (misspec): logged once per fit, the same in every outcome column: every diagnostic but nnz_or* agrees across outcome columns
- [x] diagnostics (misspec): replications per cell = 10000: min 10000, max 10000
- [x] full (correct): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 144 arm-cells; largest |difference| between the routes
- [x] (info) full (correct): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.00, middle 90% 0.89 to 1.13
- [x] full (misspec): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree: 192 arm-cells; largest |difference| between the routes
- [x] (info) full (misspec): MAE MCSE against a plain bootstrap: stored / plain bootstrap (B = 100): median 1.01, middle 90% 0.90 to 1.13
- [x] diagnostics (correct): medians, quartiles, SD, range agree: 256 rows; largest |difference| between the routes
- [x] diagnostics (correct): path truncation counts agree: 
- [x] diagnostics (misspec): medians, quartiles, SD, range agree: 256 rows; largest |difference| between the routes
- [x] diagnostics (misspec): path truncation counts agree: 

## app_p_approx (passed)

Built 2026-10-09 00:45; files: aipw_correct_p_dgp2.csv.gz, aipw_diag_correct_p_dgp2.csv.gz, aipw_misspec_p_dgp2.csv.gz, aipw_diag_misspec_p_dgp2.csv.gz.

- [x] add-on pairing (correct): oracle_hajek in aipw_diag_correct_p_dgp2.csv.gz = aipw_correct_p_dgp2.csv.gz: 360000 of 360000 rows matched, max |diff| 0
- [x] approximation bias (correct): abw_balnet0 and ols on the same draws, each cell once: 12 of 36 cells kept
- [x] add-on pairing (misspec): oracle_hajek in aipw_diag_misspec_p_dgp2.csv.gz = aipw_misspec_p_dgp2.csv.gz: 360000 of 360000 rows matched, max |diff| 0
- [x] approximation bias (misspec): abw_balnet0 and ols on the same draws, each cell once: 12 of 36 cells kept
- [x] approximation bias (correct): per-cell biases and MCSEs agree: 
- [x] approximation bias (correct): averages over cells and MCSEs agree: 
- [x] approximation bias (misspec): per-cell biases and MCSEs agree: 
- [x] approximation bias (misspec): averages over cells and MCSEs agree: 

## app_misspec_nfine (passed)

Built 2026-10-09 00:46; files: aipw_full_misspec_nfine_dgp2.csv.gz.

- [x] approximation bias (misspec): abw_balnet0 and ols on the same draws, each cell once: 87 of 87 cells kept
- [x] approximation bias (misspec): per-cell biases and MCSEs agree: 
- [x] approximation bias (misspec): averages over cells and MCSEs agree: 

## cv_overlap (passed)

Built 2026-10-12 18:57; files: r2500/cv_overlap_dgp2.csv.gz.

- [x] main: 2500 replications of every row in every cell: 6 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 0 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 2.7e-15 (IPW-balnet), 1.4e-15 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 6 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 1296 recomputed, 1296 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 3.3e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 36 recomputed, 36 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 3.3e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 72 recomputed, 72 stored, largest |difference| 8.5e-18
- [x] cv pairs: every pair in every cell found in both routes: 1296 recomputed, 1296 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 4.7e-16
- [x] cv sums: every set found in both routes: 333 recomputed, 333 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 8.9e-16
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 6 design cells: median 1.00, all within 0.96 to 1.11 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 24 sampled cells: path end RMSE median 0.97, log2_R_best median 0.99, all within 0.86 to 1.11 (B = 200)

## cv_noise (passed)

Built 2026-10-12 18:57; files: r2500/cv_noise_dgp2.csv.gz.

- [x] main: 2500 replications of every row in every cell: 7 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 0 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 1.7e-15 (IPW-balnet), 1.0e-15 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 7 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 1512 recomputed, 1512 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 2.2e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 42 recomputed, 42 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 2.6e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 84 recomputed, 84 stored, largest |difference| 9.0e-18
- [x] cv pairs: every pair in every cell found in both routes: 1512 recomputed, 1512 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 5.9e-16
- [x] cv sums: every set found in both routes: 378 recomputed, 378 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 8.9e-16
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 7 design cells: median 1.02, all within 0.92 to 1.06 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 28 sampled cells: path end RMSE median 0.98, log2_R_best median 1.01, all within 0.87 to 1.09 (B = 200)

## cv_noise_x_overlap (passed)

Built 2026-10-12 18:59; files: r2500/cv_noise_x_overlap_dgp2.csv.gz.

- [x] main: 2500 replications of every row in every cell: 28 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 0 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 2.7e-15 (IPW-balnet), 1.0e-15 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 28 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 6048 recomputed, 6048 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 4.4e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 168 recomputed, 168 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 4.4e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 336 recomputed, 336 stored, largest |difference| 8.0e-18
- [x] cv pairs: every pair in every cell found in both routes: 6048 recomputed, 6048 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 5.3e-16
- [x] cv sums: every set found in both routes: 1530 recomputed, 1530 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 1.1e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 28 design cells: median 1.00, all within 0.88 to 1.11 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 30 sampled cells: path end RMSE median 1.00, log2_R_best median 1.00, all within 0.89 to 1.18 (B = 200)

## cv_misspec (passed)

Built 2026-10-12 19:02; files: r2500/cv_misspec_dgp2.csv.gz.

- [x] main: 2500 replications of every row in every cell: 6 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 0 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 3.0e-15 (IPW-balnet), 1.0e-15 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 6 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 9072 recomputed, 9072 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 3.3e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 252 recomputed, 252 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 3.3e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 72 recomputed, 72 stored, largest |difference| 7.0e-18
- [x] cv pairs: every pair in every cell found in both routes: 9072 recomputed, 9072 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 5.6e-16
- [x] cv sums: every set found in both routes: 1080 recomputed, 1080 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 2.7e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 6 design cells: median 0.96, all within 0.94 to 1.08 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 24 sampled cells: path end RMSE median 1.00, log2_R_best median 0.99, all within 0.87 to 1.12 (B = 200)

## cv_folds (passed)

Built 2026-10-12 19:05; files: r2500/cv_folds_k10_dgp2.csv.gz, r2500/cv_folds_k3_dgp2.csv.gz, r2500/cv_folds_k5_dgp2.csv.gz.

- [x] k10: 2500 replications of every row in every cell: 4 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] k10: replications that failed (err not empty): 0 rows
- [x] k10: path_30 equals the fixed arm in every draw: largest difference 4.0e-15 (IPW-balnet), 1.0e-15 (ABW)
- [x] k10: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] k10: a lambda grid of 30 points for every cell: 4 cells in the grid file
- [x] k3: 2500 replications of every row in every cell: 4 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] k3: replications that failed (err not empty): 0 rows
- [x] k3: path_30 equals the fixed arm in every draw: largest difference 4.0e-15 (IPW-balnet), 1.0e-15 (ABW)
- [x] k3: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] k3: a lambda grid of 30 points for every cell: 4 cells in the grid file
- [x] k3: path and path end equal the first file's in every draw: 620000 rows, largest difference 0.0e+00
- [x] k5: 2500 replications of every row in every cell: 4 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] k5: replications that failed (err not empty): 0 rows
- [x] k5: path_30 equals the fixed arm in every draw: largest difference 4.0e-15 (IPW-balnet), 1.0e-15 (ABW)
- [x] k5: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] k5: a lambda grid of 30 points for every cell: 4 cells in the grid file
- [x] k5: path and path end equal the first file's in every draw: 620000 rows, largest difference 0.0e+00
- [x] cv perf: every arm-cell found in both routes: 7056 recomputed, 7056 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 3.3e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 168 recomputed, 168 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 3.3e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 96 recomputed, 96 stored, largest |difference| 9.0e-18
- [x] cv pairs: every pair in every cell found in both routes: 7728 recomputed, 7728 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 5.6e-16
- [x] cv sums: every set found in both routes: 1482 recomputed, 1482 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 2.7e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 4 design cells: median 1.02, all within 0.92 to 1.04 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 16 sampled cells: path end RMSE median 0.98, log2_R_best median 1.01, all within 0.89 to 1.10 (B = 200)

## cv_sparsity (passed)

Built 2026-10-12 19:12; files: r2500/cv_sparsity_dgp3.csv.gz.

- [x] main: 2500 replications of every row in every cell: 9 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 0 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 0.0e+00 (IPW-balnet), 0.0e+00 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 9 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 13608 recomputed, 13608 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 3.3e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 378 recomputed, 378 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 4.4e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 108 recomputed, 108 stored, largest |difference| 0.0e+00
- [x] cv pairs: every pair in every cell found in both routes: 13608 recomputed, 13608 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 6.1e-16
- [x] cv sums: every set found in both routes: 1584 recomputed, 1584 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 2.0e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 9 design cells: median 1.02, all within 0.93 to 1.07 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 30 sampled cells: path end RMSE median 0.99, log2_R_best median 0.99, all within 0.88 to 1.15 (B = 200)

## cv_sparsity_misspec (passed)

Built 2026-10-12 19:20; files: r2500/cv_sparsity_misspec_dgp3.csv.gz.

- [x] main: 2500 replications of every row in every cell: 9 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 0 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 0.0e+00 (IPW-balnet), 0.0e+00 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 9 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 13608 recomputed, 13608 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 4.4e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 378 recomputed, 378 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 4.4e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 108 recomputed, 108 stored, largest |difference| 0.0e+00
- [x] cv pairs: every pair in every cell found in both routes: 13608 recomputed, 13608 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 5.6e-16
- [x] cv sums: every set found in both routes: 1593 recomputed, 1593 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 1.8e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 9 design cells: median 0.98, all within 0.92 to 1.02 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 30 sampled cells: path end RMSE median 0.99, log2_R_best median 1.01, all within 0.85 to 1.16 (B = 200)

## cv_nbig (passed)

Built 2026-10-12 19:22; files: r2500/cv_n5000_dgp2.csv.gz, r1000/cv_n10000_dgp2.csv.gz.

- [x] n5000: 2500 replications of every row in every cell: 2 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] n5000: replications that failed (err not empty): 0 rows
- [x] n5000: path_30 equals the fixed arm in every draw: largest difference 2.0e-15 (IPW-balnet), 1.5e-15 (ABW)
- [x] n5000: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] n5000: a lambda grid of 30 points for every cell: 2 cells in the grid file
- [x] n10000: 1000 replications of every row in every cell: 2 cells, 108 estimator rows, n per cell 1000 to 1000
- [x] n10000: replications that failed (err not empty): 0 rows
- [x] n10000: path_30 equals the fixed arm in every draw: largest difference 2.0e-15 (IPW-balnet), 1.0e-15 (ABW)
- [x] n10000: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] n10000: a lambda grid of 30 points for every cell: 2 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 6048 recomputed, 6048 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 3.3e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 168 recomputed, 168 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 4.4e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 48 recomputed, 48 stored, largest |difference| 5.0e-18
- [x] cv pairs: every pair in every cell found in both routes: 6048 recomputed, 6048 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 5.8e-16
- [x] cv sums: every set found in both routes: 684 recomputed, 684 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 2.7e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 4 design cells: median 1.00, all within 0.92 to 1.05 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 16 sampled cells: path end RMSE median 1.02, log2_R_best median 1.03, all within 0.88 to 1.09 (B = 200)

## cv_alpha_x_overlap (passed)

Built 2026-10-12 19:34; files: r2500/cv_alpha_x_overlap_dgp2.csv.gz.

- [x] main: 2500 replications of every row in every cell: 20 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 0 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 6.0e-15 (IPW-balnet), 3.0e-15 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 20 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 30240 recomputed, 30240 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 3.3e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 840 recomputed, 840 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 4.4e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 240 recomputed, 240 stored, largest |difference| 3.5e-18
- [x] cv pairs: every pair in every cell found in both routes: 30240 recomputed, 30240 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 6.2e-16
- [x] cv sums: every set found in both routes: 3501 recomputed, 3501 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 3.6e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 20 design cells: median 1.00, all within 0.89 to 1.06 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 30 sampled cells: path end RMSE median 1.00, log2_R_best median 1.02, all within 0.89 to 1.09 (B = 200)

## cv_n (passed)

Built 2026-10-12 19:46; files: r2500/cv_n_dgp2.csv.gz.

- [x] main: 2500 replications of every row in every cell: 16 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 0 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 0.0e+00 (IPW-balnet), 0.0e+00 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 16 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 24192 recomputed, 24192 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 3.3e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 672 recomputed, 672 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 4.4e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 192 recomputed, 192 stored, largest |difference| 0.0e+00
- [x] cv pairs: every pair in every cell found in both routes: 24192 recomputed, 24192 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 5.6e-16
- [x] cv sums: every set found in both routes: 2808 recomputed, 2808 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 2.2e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 16 design cells: median 1.03, all within 0.96 to 1.12 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 30 sampled cells: path end RMSE median 1.00, log2_R_best median 1.01, all within 0.88 to 1.19 (B = 200)

## cv_p (passed)

Built 2026-10-12 19:58; files: r2500/cv_p_dgp2.csv.gz.

- [x] main: 2500 replications of every row in every cell: 16 cells, 108 estimator rows, n per cell 2500 to 2500
- [x] main: replications that failed (err not empty): 108 rows
- [x] main: path_30 equals the fixed arm in every draw: largest difference 2.6e-15 (IPW-balnet), 1.0e-15 (ABW)
- [x] main: chosen lambdas the same in every error column: largest spread 0.0e+00
- [x] main: a lambda grid of 30 points for every cell: 16 cells in the grid file
- [x] cv perf: every arm-cell found in both routes: 24192 recomputed, 24192 stored
- [x] cv perf: bias, EmpSE, RMSE, share and MCSEs agree: largest |difference| between the routes 4.4e-16
- [x] cv perf: replication counts agree: 0 arm-cells differ
- [x] cv cells: every cell found in both routes: 672 recomputed, 672 stored
- [x] cv cells: RMSEs, MCSEs, ratios and log2 ratios agree: largest |difference| between the routes 3.9e-16
- [x] cv cells: path minimum, best selector and the three rules agree: 0 cells where the routes choose differently
- [x] cv picks: median chosen lambda agrees: 192 recomputed, 192 stored, largest |difference| 5.5e-18
- [x] cv pairs: every pair in every cell found in both routes: 24192 recomputed, 24192 stored
- [x] cv pairs: log R_c and its paired MCSE agree: largest |difference| between the routes 6.1e-16
- [x] cv sums: every set found in both routes: 2772 recomputed, 2772 stored
- [x] cv sums: sums of logs and their MCSEs agree: largest |difference| between the routes 3.6e-15
- [x] (info) cv bootstrap/delta MCSE (sum over a design cell's columns): bootstrap SD / delta MCSE over 16 design cells: median 1.00, all within 0.92 to 1.09 (B = 200)
- [x] (info) cv bootstrap/delta MCSE (path end RMSE, log2_R_best): bootstrap SD / delta MCSE over 30 sampled cells: path end RMSE median 0.99, log2_R_best median 0.98, all within 0.89 to 1.18 (B = 200)

## wc_compare (passed)

Built 2026-10-12 22:11; files: wc_noise_s1_c1.rds, wc_noise_s1_c2.rds, wc_noise_s1_c4.rds, wc_noise_s10_c1.rds, wc_noise_s10_c2.rds, wc_noise_s10_c4.rds, wc_noise_s30_c1.rds, wc_noise_s30_c2.rds, wc_noise_s30_c4.rds, wc_noise_s100_c1.rds, wc_noise_s100_c2.rds, wc_noise_s100_c4.rds.

- [x] s1_c1: 1000 replications of every arm, none failed: 1000 run, 0 failed, 106 arms per outcome
- [x] s1_c1: path endpoint reached by 95% of replications: endpoint at grid point 100 of 100
- [x] s1_c1: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s1_c2: 1000 replications of every arm, none failed: 1000 run, 0 failed, 106 arms per outcome
- [x] s1_c2: path endpoint reached by 95% of replications: endpoint at grid point 100 of 100
- [x] s1_c2: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s1_c4: 1000 replications of every arm, none failed: 1000 run, 0 failed, 57 arms per outcome
- [x] s1_c4: path endpoint reached by 95% of replications: endpoint at grid point 51 of 100
- [x] s1_c4: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s10_c1: 1000 replications of every arm, none failed: 1000 run, 0 failed, 106 arms per outcome
- [x] s10_c1: path endpoint reached by 95% of replications: endpoint at grid point 100 of 100
- [x] s10_c1: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s10_c2: 1000 replications of every arm, none failed: 1000 run, 0 failed, 106 arms per outcome
- [x] s10_c2: path endpoint reached by 95% of replications: endpoint at grid point 100 of 100
- [x] s10_c2: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s10_c4: 1000 replications of every arm, none failed: 1000 run, 0 failed, 57 arms per outcome
- [x] s10_c4: path endpoint reached by 95% of replications: endpoint at grid point 51 of 100
- [x] s10_c4: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s30_c1: 1000 replications of every arm, none failed: 1000 run, 0 failed, 106 arms per outcome
- [x] s30_c1: path endpoint reached by 95% of replications: endpoint at grid point 100 of 100
- [x] s30_c1: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s30_c2: 1000 replications of every arm, none failed: 1000 run, 0 failed, 106 arms per outcome
- [x] s30_c2: path endpoint reached by 95% of replications: endpoint at grid point 100 of 100
- [x] s30_c2: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s30_c4: 1000 replications of every arm, none failed: 1000 run, 0 failed, 57 arms per outcome
- [x] s30_c4: path endpoint reached by 95% of replications: endpoint at grid point 51 of 100
- [x] s30_c4: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s100_c1: 1000 replications of every arm, none failed: 1000 run, 0 failed, 106 arms per outcome
- [x] s100_c1: path endpoint reached by 95% of replications: endpoint at grid point 100 of 100
- [x] s100_c1: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s100_c2: 1000 replications of every arm, none failed: 1000 run, 0 failed, 106 arms per outcome
- [x] s100_c2: path endpoint reached by 95% of replications: endpoint at grid point 100 of 100
- [x] s100_c2: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] s100_c4: 1000 replications of every arm, none failed: 1000 run, 0 failed, 57 arms per outcome
- [x] s100_c4: path endpoint reached by 95% of replications: endpoint at grid point 51 of 100
- [x] s100_c4: the five shared selectors present: cv.bloss, cv.smd, cv.inf, boot.smd, boot.inf
- [x] DGP 2 twin: the stored summary passed its own checks: cv_misspec built 2026-10-12 19:02
- [x] DGP 2 twin: every matched overlap level present: 0.25, 0.5, 1
- [x] every Wong-Chan point matched to a DGP 2 bracket: 15 matched pairs
- [x] wc points: every Wong-Chan point found in both routes: 24 recomputed, 24 stored
- [x] wc points: log ratios, their paired MCSEs and the endpoint RMSE agree: largest |difference| between the routes 2.8e-16
- [x] wc points: endpoint, best fixed lambda, best selector and the two rules agree: 0 points where the routes choose differently
- [x] wc dgp2: every DGP 2 point found in both routes: 21 recomputed, 21 stored
- [x] wc dgp2: log ratios from the arms' RMSEs equal the paired ones: largest |difference| 0.0e+00
- [x] wc dgp2: best fixed lambda, best selector and the two rules agree: 0 points where the routes choose differently
- [x] (info) wc dgp2: best selector as in the twin summary's cells table: 0 points differ
- [x] wc dgp2: the twin summary is the one the entry was built from: cv_misspec md5 8ac00712
- [x] wc twins: the bracketing cells and both agreement flags agree: 0 of 15 pairs differ
- [x] (info) wc draft: four points as the 9 October draft reported them: largest |difference| to six decimals 4.9e-07
- [x] (info) wc bootstrap/delta MCSE (best selector's log ratio): bootstrap SD / delta MCSE at 3 points: 1.02, 1.01, 1.02 (B = 200)

## Across sections

- [x] (info) path ends, sparsity grid (s = 4, nu = 1, iid) against the p/n grid: correct score, p/n 0.5 and 1 (4 cells): largest |difference| 0.42% (treated or control arm), 0.32% (their sum); largest |z| 1.78 (independent draws)

