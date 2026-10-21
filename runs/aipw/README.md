# runs/aipw

One script per design axis. Each `run_arm()` call writes one file to
`results/aipw/` and is skipped if that file exists, so sourcing a script
only runs what is missing. Shared code is in `setup.R`.

| Script | Results files |
|---|---|
| `n_axis.R` | `aipw_{correct,misspec}_n`, `aipw_add_misspec_n`, `aipw_tol_{correct,misspec}_n`, `aipw_correct_n400k`, `aipw_add_correct_n400k`, `aipw_misspec_n_dgp2n1mill`, `aipw_diag_{correct,misspec}_nfine`, `aipw_full_misspec_nfine`, `aipw_{correct,misspec}_n250` |
| `overlap_axis.R` | `aipw_{correct,misspec}_overlap7`, `_overlap8`, `_n_by_overlap`, `_overlap` (4 levels), `aipw_diag_correct_overlap7`, `aipw_full_misspec_overlap7` |
| `p_axis.R` | `aipw_{correct,misspec}_p`, `aipw_diag_{correct,misspec}_p` |
| `sparsity_axis.R` | `aipw_correct_sparsity_hd_dgp2`, `aipw_misspec_tanobs_sparsity_hd_dgp3`, `aipw_{correct,misspec_tanobs}_sparsity_hd_equi_dgp3` |
| `strength_axis.R` | `aipw_{correct,misspec}_strength` |
| `treat_prop_axis.R` | `aipw_{correct,misspec}_treat_prop` |

All `_dgp2` unless marked `_dgp3`.

Run one arm only: `ARMS=tol Rscript runs/aipw/n_axis.R`.
Smoke run: `N_SIM=2 OUT_DIR=<temp folder> Rscript runs/aipw/<script>`.

## Old script names

The code now lives in:

- `aipw_correct_misspec_n_dgp2.R`, `aipw_tol_dgp2.R`, `aipw_largeN.R`,
  `aipw_correct_misspec_n1mill_dgp2.R`, `aipw_correct_misspec_n2_dgp2.R`,
  `aipw_correct_misspec_n250to200k.R`, `aipw_n250_dgp2.R` → `n_axis.R`
- `aipw_overlap.R`, `aipw_overlap10k_dgp2.R`, `aipw_n_by_overlap_dgp2.R`,
  `aipw_misspec_overlap.R` → `overlap_axis.R`
- `aipw_addon_dgp2.R` (one-off) → its arms are in `n_axis.R` and `overlap_axis.R`
- `aipw_paxis_dgp2.R`, `aipw_diag_p_dgp2.R` → `p_axis.R`
- `aipw_sparsity_hd.R`, `aipw_misspec_tanobs_sparsity_hd_dgp3.R`,
  `equicorrelation.R` → `sparsity_axis.R`
- `aipw_strength.R` → `strength_axis.R`
- `aipw_treat_prop.R` → `treat_prop_axis.R`
