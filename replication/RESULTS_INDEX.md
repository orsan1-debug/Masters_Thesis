# Results index

Every table and figure of the thesis, the stored summary it is drawn from
(`output/thesis_summaries/<id>.rds`), the run that reproduces the raw results and
the stored files it reproduces. Stored Study 1 files are in `results/aipw/`, Study 2
files in `results/cv/our dgp/r2500/` (`r1000/` for the 1000-replication run) and the
Wong and Chan cells in `results/cv/wz_replication/`. A run is one grid and writes one file,
`results/study1/<grid>.csv.gz` or `results/study2/<grid>.csv.gz`, with every
estimator. The stored files split some grids by specification, fold count or
estimator set over separate runs that shared their draws.

## Study 1 (`runs/study1.R`, base seed 202)

| Run | Grid | Stored files it reproduces | Stored summaries | Thesis displays |
|---|---|---|---|---|
| `sample_size` | n = 500 to 400 000, c = 1, 10 000 replications. Misspecified to n = 50 000 | Correct: `aipw_correct_n400k_dgp2`, `aipw_correct_n_dgp2`, `aipw_add_correct_n400k_dgp2`. Misspecified: `aipw_misspec_n_dgp2`, `aipw_add_misspec_n_dgp2` | `correct_n`, `misspec_n`, `app_misspec_n`, `n_both` | Table 4, Appendix Tables 1, 2, 4 and 5 |
| `sample_size_path` | n = 250 to 200 000 in 29 steps and n = 1 000 000, misspecified | `aipw_full_misspec_nfine_dgp2`, `aipw_misspec_n_dgp2n1mill` (stored as ten batches of 1000) | `misspec_npath` | Figures 1–2, Appendix Table 3 |
| `overlap` | c = 0.05 to 2 at n = 1000 and n = 10 000 | Correct: `aipw_correct_overlap7_dgp2`, `aipw_diag_correct_overlap7_dgp2`, `aipw_correct_overlap8_dgp2`. Misspecified: `aipw_misspec_overlap7_dgp2`, `aipw_full_misspec_overlap7_dgp2`, `aipw_misspec_overlap8_dgp2` | `overlap`, `app_overlap_full` | Tables 5–7, Figure 3, Appendix Tables 7–8, Appendix Figure 2 |
| `n_by_overlap` | n = 100 to 5000 by c = 0.01 to 1, 1000 replications | `aipw_correct_n_by_overlap_dgp2`, `aipw_misspec_n_by_overlap_dgp2` | `nc`, `app_nc` | Appendix Figure 3, Appendix Tables 9–10 |
| `dimension` | p/n = 0.05 to 4 at n = 500 and 1000 | Correct: `aipw_correct_p_dgp2`, `aipw_diag_correct_p_dgp2`. Misspecified: `aipw_misspec_p_dgp2`, `aipw_diag_misspec_p_dgp2` | `p_axis`, `app_p_axis` | Table 8, Figure 4, Appendix Tables 11–14, Appendix Figure 4 |
| `outcome_signal` | kappa = 0.02 to 1 at n = 1000 | `aipw_correct_strength_dgp2`, `aipw_misspec_strength_dgp2` | `signal`, `app_signal` | Table 9, Appendix Tables 21–22, Appendix Figure 5 |
| `denser` | s = 4, 16, 64 by nu, p/n, n. Independent, AR(1) and equicorrelated, 1000 replications | Correct: `aipw_correct_sparsity_hd_dgp2`, `aipw_correct_sparsity_hd_equi_dgp3`. Misspecified: `aipw_misspec_tanobs_sparsity_hd_dgp3`, `aipw_misspec_tanobs_sparsity_hd_equi_dgp3` | `sparsity`, `app_sparsity`, `app_sparsity_equi`, `app_sparsity_tan` | Figure 5, Appendix Tables 15–19 |

Figure 6, Appendix Figures 2 and 5 and the coverage text (Results 4.5) read
`coverage_empse`, built by the qmd chunk `cov-data` from the stored files of the
`sample_size`, `sample_size_path` (n = 1 000 000), `overlap` and `outcome_signal`
runs. Appendix Figure 1, Appendix Table 6 and Appendix Table 20 come from the design
alone. The qmd computes them in their own chunks, and `share_outside()` and
`common_factor_share()` in `R/dgp.R` give the same two tables (`tests/test_dgp.R`).

The summaries `app_correct_n`, `app_misspec_nfine`, `app_overlap`, `app_p_approx`,
`app_ipw_800k`, `treat_share`, `tol_n` and `cv_folds` are stored but no display
reads them. No display reads `aipw_correct_n_dgp2`, the `aipw_diag_` files or
`aipw_full_misspec_overlap7_dgp2`.

## Study 2 (`runs/study2.R`)

| Run | Grid | Seed | Stored files it reproduces | Stored summary | Thesis displays |
|---|---|---|---|---|---|
| `overlap` | c = 1 to 0.05 | 301 | `cv_overlap_dgp2` | `cv_overlap` | Table 10, Appendix Table 24 |
| `noise` | kappa = 1 to 0 | 302 | `cv_noise_dgp2` | `cv_noise` | as above |
| `noise_by_overlap` | kappa by c | 303 | `cv_noise_x_overlap_dgp2` | `cv_noise_x_overlap` | Table 10, Figures 7–8, Appendix Table 24, Appendix Figure 8 |
| `penalty_by_overlap` | alpha = 0 to 1 by c | 305 | `cv_alpha_x_overlap_dgp2` | `cv_alpha_x_overlap` | Table 10, Appendix Tables 24–25, Appendix Figure 7 |
| `sample_size` | n = 250 to 2000 by c by alpha | 306 | `cv_n_dgp2` | `cv_n` | Table 10, Appendix Table 24, Appendix Figure 7 |
| `dimension` | p = 50 to 1000 by c by alpha | 307 | `cv_p_dgp2` | `cv_p` | as above |
| `misspecification` | c = 1 to 0.05, misspecified | 308 | `cv_misspec_dgp2` | `cv_misspec` | Table 10, Appendix Table 24, Appendix Figure 8, Appendix H |
| `folds` | c = 1 to 0.15 on 10, 5 and 3 folds | 309 | 10 folds: `cv_folds_k10_dgp2`. 5: `cv_folds_k5_dgp2`. 3: `cv_folds_k3_dgp2` | `cv_folds` | none, run but not reported |
| `large_n_5000` | n = 5000 by c, 2500 replications | 310 | `cv_n5000_dgp2` | `cv_nbig` | Table 10, Appendix Table 24, Appendix Figure 7 |
| `large_n_10000` | n = 10 000 by c, 1000 replications | 311 | `r1000/cv_n10000_dgp2` | `cv_nbig` | as above |
| `denser` | s = 4, 16, 64 by covariance, correct and under Tan's transformation | 312 | Correct: `cv_sparsity_dgp3`. Misspecified: `cv_sparsity_misspec_dgp3` | `cv_sparsity`, `cv_sparsity_misspec` | Table 10, Appendix Table 24, Appendix Figure 8 |

Each Study 2 results file, stored or rerun, has a `_lamgrid.csv` beside it with the
fixed lambda grid of every cell. Every Study 1 and Study 2 results file has a
`_session.txt` with the session that produced it. The Wong and Chan files store their
R, balnet and glmnet versions inside.

Table 3 gives the large n grid 1000 replications. The n = 5000 run has 2500, as
stored, and the thesis's Study 2 numbers at n = 5000 come from it.

## Wong and Chan check (`runs/wong_chan.R`, master seed 20260903)

Sixteen cells, noise sd 1, 10, 30, 100 by logit multiplier 1, 2, 3, 4, one
`wc_noise_s<sigma>_c<c>.rds` each. The thesis matches the cells with multipliers
1, 2 and 4 to the `misspecification` grid above (stored summary `wc_compare`,
Appendix H, Appendix Table 23, Appendix Figure 6).

## Design tables

Tables 1–3 are typed in the qmd from the design. The grids in `runs/study1.R` and
`runs/study2.R` are the same grids in code. Table 2 also lists the treated
proportion and path endpoint grids, which were run but are not reported and are not
in `runs/study1.R`. `runs/study2.R` also has the folds grid, which Table 3 does not
list.

## Displays and qmd cells

Each display of the thesis, the cell of `thesis/AAAAA Working Thesis.qmd` that draws
it and the stored summaries that cell reads. The ten Study 2 summaries are the ten in
the Study 2 table above without `cv_folds`.

| Display | qmd cell | Reads |
|---|---|---|
| Table 1 | `tbl-spec` | typed |
| Table 2 | `tbl-grids` | typed |
| Table 3 | `tbl-grids-s2` | typed |
| Table 4 | `tbl-n-adv` | `n_both` |
| Figure 1 | `fig-n-refit` | `misspec_npath` |
| Figure 2 | `fig-n-misspec` | `misspec_npath` |
| Table 5 | `tbl-ov-adv` | `overlap` |
| Table 6 | `tbl-ov-adv-hi` | `overlap` |
| Table 7 | `tbl-ov-adv-lo` | `overlap` |
| Figure 3 | `fig-ov-adv` | `overlap` |
| Table 8 | `tbl-pn-adv` | `p_axis` |
| Figure 4 | `fig-pn-rmse` | `p_axis`, `app_p_axis` |
| Figure 5 | `fig-sp-change` | `sparsity` |
| Table 9 | `tbl-sig-adv` | `signal` |
| Figure 6 | `fig-n-coverage` | `coverage_empse` (built in `cov-data`) |
| Table 10 | `tbl-s2-counts` | the ten Study 2 summaries |
| Figure 7 | `fig-s2-path` | `cv_noise_x_overlap` |
| Figure 8 | `fig-s2-lead` | `cv_noise_x_overlap` |
| Table 11 | `tbl-concl` | `n_both`, `overlap`, `p_axis`, `sparsity`, `signal` |
| Appendix Table 1 | `atbl-q1-adv` | `correct_n` |
| Appendix Table 2 | `atbl-app-n-diag-m-lin` | `app_misspec_n` |
| Appendix Table 3 | `atbl-n-path` | `misspec_npath` |
| Appendix Table 4 | `atbl-app-avg-n-m` | `misspec_n` |
| Appendix Table 5 | `atbl-app-n-ref` | `n_both` |
| Appendix Figure 1 | `afig-ov-levels` | design, computed in the cell |
| Appendix Table 6 | `atbl-ov-levels` | design, computed in the cell |
| Appendix Table 7 | `atbl-app-ov-full-c` | `app_overlap_full` |
| Appendix Table 8 | `atbl-app-ov-full-m` | `app_overlap_full` |
| Appendix Figure 2 | `afig-ov-coverage` | `coverage_empse` |
| Appendix Figure 3 | `afig-app-nc` | `nc` |
| Appendix Table 9 | `atbl-app-nc-full-m` | `app_nc` |
| Appendix Table 10 | `atbl-app-nc-trunc` | `nc` |
| Appendix Figure 4 | `afig-app-pn-bias` | `p_axis`, `app_p_axis` |
| Appendix Table 11 | `atbl-app-pn-full-c` | `app_p_axis` |
| Appendix Table 12 | `atbl-app-pn-full-m` | `app_p_axis` |
| Appendix Table 13 | `atbl-app-pn-diag-c-lin` | `app_p_axis` |
| Appendix Table 14 | `atbl-app-pn-diag-m-lin` | `app_p_axis` |
| Appendix Table 15 | `atbl-app-sp-cells` | `app_sparsity`, `app_sparsity_equi` |
| Appendix Table 16 | `atbl-app-sp-bound` | `sparsity` |
| Appendix Table 17 | `atbl-app-sp-exact` | `sparsity` |
| Appendix Table 18 | `atbl-app-sp-trunc` | `sparsity` |
| Appendix Table 19 | `atbl-app-sp-tan` | `app_sparsity_tan` |
| Appendix Table 20 | `atbl-app-sp-share` | design, computed in the cell |
| Appendix Table 21 | `atbl-app-sig-full-c` | `app_signal` |
| Appendix Table 22 | `atbl-app-sig-diag-c-lin` | `app_signal` |
| Appendix Figure 5 | `afig-sig-coverage` | `coverage_empse` |
| Appendix Figure 6 | `afig-app-wc` | `wc_compare` |
| Appendix Table 23 | `atbl-app-wc` | `wc_compare` |
| Appendix Table 24 | `atbl-s2-selectors` | the ten Study 2 summaries |
| Appendix Table 25 | `atbl-s2-alpha` | `cv_alpha_x_overlap` |
| Appendix Figure 7 | `afig-s2-design` | `cv_alpha_x_overlap`, `cv_n`, `cv_nbig`, `cv_p` |
| Appendix Figure 8 | `afig-s2-spec` | `cv_noise_x_overlap`, `cv_misspec`, `cv_sparsity`, `cv_sparsity_misspec` |
