# Pending consolidation (script-local helpers)

Twenty-two names are defined with different bodies in more than one script and
have no `R/` version. They are per-script helpers (same name, different job per
script), so they were **not merged** in Phase 3. Listed for a later decision.
Paths are current; archived scripts are marked. Diffs: `notes/dropped_variants.md`.

## `arity_ok`  (3 files, 2 distinct bodies)

- `runs/ipw/cvfix.R`
- `tests/balnet_cvfix_gate.R`
- `archive/tests/balnet_fix.R` (archived)

## `att`  (5 files, 2 distinct bodies)

- `runs/cv/wz_replication/run_wc_basis.R`
- `runs/cv/wz_replication/run_wc_noise.R`
- `runs/cv/wz_replication/run_wc_overlap.R`
- `runs/cv/wz_replication/run_wc_overlap_n.R`
- `runs/cv/wz_replication/run_wc_path.R`

## `axis_panels`  (2 files, 2 distinct bodies)

- `analysis/cv/basic_dgp/plot_grids.R`
- `archive/extension_cv/scripts/plot_rmse_grids_for_qmd.R` (archived)

## `bind`  (17 files, 3 distinct bodies)

- `analysis/cv/wz_replication/plot_wc_att.R`
- `analysis/cv/wz_replication/wc_findings.qmd`
- `analysis/cv/wz_replication/wc_grids.R`
- `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`
- `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`
- `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`
- `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- `archive/exploring_cv/balnet_mini_experiment_edited.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_machine_2.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_windows.R` (archived)

## `dgp_gen`  (13 files, 11 distinct bodies)

- `runs/ipw/correct_misspec_dgp2.R`
- `runs/ipw/correctspec_dgp2_n800k.R`
- `runs/ipw/cvfix.R`
- `runs/ipw/DimensionDGP2.R`
- `runs/ipw/DimensionDGP2N5K.R`
- `runs/ipw/overlap_misspec_dgp2.R`
- `runs/ipw/PoverMDGP2.R`
- `runs/ipw/SNR.R`
- `runs/ipw/Sparsity.R`
- `runs/ipw/sparsity_2.R`
- `runs/ipw/treat_prop.R`
- `runs/ipw/weak_instruments.R`
- `archive/tests/balnet_fix.R` (archived)

## `fit`  (28 files, 10 distinct bodies)

- `R/estimators_cv.R`
- `runs/cv/basic_dgp/run_alpha_sweep_bad_overlap__alpha_bad.R`
- `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`
- `runs/cv/basic_dgp/run_bad_overlap_noise10_1000reps__tune4_19.R`
- `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`
- `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`
- `runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R`
- `runs/cv/basic_dgp/run_noise_axis_moderate_bad__snr_mb.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- `runs/cv/basic_dgp/run_overlap_axis__ov_s1.R`
- `runs/cv/basic_dgp/run_overlap_axis_1000reps__ov_1k.R`
- `runs/cv/basic_dgp/run_overlap_axis_maxit1e5__ov_s1_m5.R`
- `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`
- `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- `runs/cv/wz_replication/run_wc_basis.R`
- `runs/cv/wz_replication/run_wc_noise.R`
- `runs/cv/wz_replication/run_wc_overlap.R`
- `runs/cv/wz_replication/run_wc_overlap_n.R`
- `runs/cv/wz_replication/run_wc_path.R`
- `tests/cv_smoke_test.R`
- `archive/exploring_cv/balnet_mini_experiment_edited.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_machine_2.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_windows.R` (archived)

## `fmt`  (3 files, 3 distinct bodies)

- `analysis/cv/basic_dgp/cv_summary.qmd`
- `analysis/cv/wz_replication/wc_grids.R`
- `archive/extension_cv/cv_summary.qmd` (archived)

## `gen_cell`  (20 files, 12 distinct bodies)

- `runs/cv/basic_dgp/run_alpha_sweep_bad_overlap__alpha_bad.R`
- `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`
- `runs/cv/basic_dgp/run_bad_overlap_noise10_1000reps__tune4_19.R`
- `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`
- `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`
- `runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R`
- `runs/cv/basic_dgp/run_noise_axis_moderate_bad__snr_mb.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- `runs/cv/basic_dgp/run_overlap_axis__ov_s1.R`
- `runs/cv/basic_dgp/run_overlap_axis_1000reps__ov_1k.R`
- `runs/cv/basic_dgp/run_overlap_axis_maxit1e5__ov_s1_m5.R`
- `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`
- `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- `archive/exploring_cv/balnet_mini_experiment_edited.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_machine_2.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_windows.R` (archived)

## `generate_data`  (2 files, 2 distinct bodies)

- `tests/dgp2_validity.R`
- `archive/tests/dgp2_validity_t1.R` (archived)

## `grid_png`  (2 files, 2 distinct bodies)

- `analysis/cv/basic_dgp/plot_grids.R`
- `archive/extension_cv/scripts/plot_rmse_grids_for_qmd.R` (archived)

## `num`  (3 files, 2 distinct bodies)

- `analysis/ipw/correctspec.qmd`
- `runs/ipw/cvfix.R`
- `archive/tests/balnet_fix.R` (archived)

## `one_rep`  (27 files, 16 distinct bodies)

- `runs/cv/basic_dgp/run_alpha_sweep_bad_overlap__alpha_bad.R`
- `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`
- `runs/cv/basic_dgp/run_bad_overlap_noise10_1000reps__tune4_19.R`
- `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`
- `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`
- `runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R`
- `runs/cv/basic_dgp/run_noise_axis_moderate_bad__snr_mb.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- `runs/cv/basic_dgp/run_overlap_axis__ov_s1.R`
- `runs/cv/basic_dgp/run_overlap_axis_1000reps__ov_1k.R`
- `runs/cv/basic_dgp/run_overlap_axis_maxit1e5__ov_s1_m5.R`
- `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`
- `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- `runs/cv/wz_replication/run_wc_basis.R`
- `runs/cv/wz_replication/run_wc_noise.R`
- `runs/cv/wz_replication/run_wc_overlap.R`
- `runs/cv/wz_replication/run_wc_overlap_n.R`
- `runs/cv/wz_replication/run_wc_path.R`
- `tests/cv_smoke_test.R`
- `archive/exploring_cv/balnet_mini_experiment_edited.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_machine_2.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_windows.R` (archived)

## `panel_rmse`  (3 files, 2 distinct bodies)

- `analysis/cv/basic_dgp/plot_grids.R`
- `runs/cv/basic_dgp/make_cv_picks.R`
- `archive/extension_cv/scripts/plot_rmse_grids_for_qmd.R` (archived)

## `patch_cv_balnet`  (3 files, 2 distinct bodies)

- `runs/ipw/cvfix.R`
- `tests/balnet_cvfix_gate.R`
- `archive/tests/balnet_fix.R` (archived)

## `rewrite`  (3 files, 2 distinct bodies)

- `runs/ipw/cvfix.R`
- `tests/balnet_cvfix_gate.R`
- `archive/tests/balnet_fix.R` (archived)

## `rmse`  (19 files, 3 distinct bodies)

- `analysis/cv/wz_replication/plot_wc_att.R`
- `analysis/cv/wz_replication/wc_findings.qmd`
- `analysis/cv/wz_replication/wc_grids.R`
- `runs/cv/basic_dgp/run_alpha_sweep_bad_overlap__alpha_bad.R`
- `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`
- `runs/cv/basic_dgp/run_confounder_spread_x_penalty__spread.R`
- `runs/cv/basic_dgp/run_high_p_x_alpha_x_overlap__dimhi.R`
- `runs/cv/basic_dgp/run_n_x_overlap_enet__n_enet.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- `runs/cv/basic_dgp/run_p_x_alpha_x_overlap__dima.R`
- `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- `archive/exploring_cv/balnet_mini_experiment_edited.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_machine_2.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_windows.R` (archived)

## `scan`  (3 files, 2 distinct bodies)

- `runs/ipw/cvfix.R`
- `tests/balnet_cvfix_gate.R`
- `archive/tests/balnet_fix.R` (archived)

## `selftest_cv_fix`  (2 files, 2 distinct bodies)

- `runs/ipw/cvfix.R`
- `archive/tests/balnet_fix.R` (archived)

## `summarise_file`  (2 files, 2 distinct bodies)

- `analysis/cv/basic_dgp/prep_cv_summary.R`
- `archive/extension_cv/scripts/summarise_runs_to_csv.R` (archived)

## `target`  (3 files, 2 distinct bodies)

- `runs/ipw/cvfix.R`
- `tests/balnet_cvfix_gate.R`
- `archive/tests/balnet_fix.R` (archived)

## `title_of`  (11 files, 9 distinct bodies)

- `runs/cv/basic_dgp/run_alpha_x_overlap_x_density_50reps__tunea.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_diluted__tune.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_density_fixed__tune2.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_enet__snr_enet.R`
- `runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R`
- `runs/cv/basic_dgp/run_p_x_overlap_x_density__tunep2.R`
- `runs/cv/basic_dgp/run_ridge_enet_lowered_floor__tunea4.R`
- `runs/cv/basic_dgp/run_ridge_enet_x_overlap_x_density__tunea3.R`
- `archive/exploring_cv/balnet_mini_experiment_edited.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_machine_2.R` (archived)
- `archive/exploring_cv/balnet_mini_experiment_windows.R` (archived)

## `walk`  (3 files, 2 distinct bodies)

- `runs/ipw/cvfix.R`
- `tests/balnet_cvfix_gate.R`
- `archive/tests/balnet_fix.R` (archived)

