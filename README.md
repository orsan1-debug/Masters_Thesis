# Masters thesis: balancing weights and cross-validated tuning

Simulation study of covariate-balancing weights (balnet) against
propensity-score weighting, and of cross-validation rules for choosing the
balance tolerance. Two components share one codebase:

- `ipw`: balnet versus glmnet ATE estimators on DGP1 (legacy) and DGP2.
- `cv`: tuning-rule experiments; `wz_replication` (Wong and Chan design,
  Wang and Zubizarreta 2020), `basic_dgp` (high-dimensional treated-mean
  DGP) and `thesis_dgp` (reserved).

Everything is run from the project root (open the `.Rproj`); paths are
built with `here::here()`.

## Layout

```
R/            shared functions only, no side effects (packages.R is the exception)
  packages.R      library() calls and the select/filter aliases; source first
  dgp.R           dgp1(), dgp2(), dgp_wc(), dgp_wc_overlap(), gen_data(), run_par()
  estimators_ipw.R  ATE estimators and estimate_all()
  estimators_cv.R   sbw fits and the Wang and Zubizarreta criterion
  simulate.R      seed streams and the checkpointed grid driver
  registry.R      run ledger (registry.csv)
  summarise.R     performance measures with MCSE, diagnostics summaries
  plots.R         figure and table helpers (needs summarise.R sourced first)
runs/         thin run scripts: runs/ipw/, runs/cv/{wz_replication,basic_dgp,thesis_dgp}/
results/      raw simulation output, read-only, git-ignored (same subfolders)
analysis/     .qmd and helper scripts that read results/ and write output/
output/       generated figures, tables, summaries and rendered documents; disposable
thesis/       chapters
tests/        smoke and validity checks
notes/        manifest, protocol, phase reports
readings/     PDFs, git-ignored
archive/      superseded scripts and results; nothing here is read by live code
```

`notes/restructure_manifest.csv` records where every file came from and
`notes/rename_map.csv` every rename.

## Running a batch

From an R session at the project root:

```r
source("runs/ipw/correctspecDGP1.R")
```

Every run script sources `R/packages.R` first, then the R/ files it needs.
The replication count and the output folder come from the environment so a
smoke run never touches `results/`:

```bash
N_SIM=2 OUT_DIR=/tmp/smoke Rscript runs/ipw/correctspecDGP1.R
N_SIM=2 OUT_DIR=/tmp/smoke Rscript runs/cv/basic_dgp/run_noise_x_overlap_lasso__tune4.R
```

`N_SIM` defaults to each script's production count; `OUT_DIR` defaults to
the script's folder under `results/`. ipw batches go through
`run_batch()`, which writes `<batch>.csv.gz`, `<batch>_cvcurves.rds` and
`<batch>_session.txt` plus per-cell checkpoints that let an interrupted run
resume. Completed batches are recorded in `registry.csv`.

## Rendering

One document:

```bash
quarto render analysis/ipw/overlap.qmd
```

All documents (`_quarto.yml` lists `analysis/**` and `thesis/**`):

```bash
quarto render
```

Rendered files go to `output/rendered/`. Figures included by the cv
documents are produced by `analysis/cv/basic_dgp/plot_grids.R` and
`analysis/cv/wz_replication/plot_wc_grids.R`; summary tables by
`analysis/cv/basic_dgp/prep_cv_summary.R`.
