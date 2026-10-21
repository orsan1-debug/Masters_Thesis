# Replication package

Code to reproduce every number, table and figure in *Balancing weights against
maximum likelihood propensity scores: a simulation study of balnet across sample
size, overlap, dimension and penalty choice* (ETC5860 Masters thesis, Monash
University, 2026).

## Quick start

Open the `.Rproj`, then `Rscript main.R` from the project root. It installs the
packages, runs the tests and renders the thesis from the stored summaries. About 10
minutes, no raw results needed. The options at its top set the replications for reruns.

## Data availability

All data are simulated by the code in this package. The design is in
`replication/R/dgp.R` and the seeds in `replication/RESULTS_INDEX.md`. The stored
summaries the thesis reads are included (`output/thesis_summaries/`, 19 MB). The raw
results (`results/`, 7 GB) are not. They are free from the author on request.

## Layout

```
README.md                 this file, also as README.pdf
LICENSE.txt               the licences
main.R                    the master file. reps and steps set at its top
_quarto.yml, .here        the Quarto project and the root for here::here()
*.Rproj                   the RStudio project. open this first
.gitattributes            keeps every file byte for byte, so the manifest's md5 sums hold
.gitignore                keeps results/, session files and superseded code out
R/, runs/, handoff/       the code that ran, kept as the record
output/thesis_summaries/  the stored summaries the thesis reads
output/rendered/thesis/   the submitted PDF
thesis/                   the qmd the PDF was rendered from, and references.bib
pkg/                      balnet 0.0.4, the version that ran
replication/
  README.md               points to this README
  RESULTS_INDEX.md        every table and figure -> stored summary -> results file -> run
  VERIFICATION.md         every table and figure, its sources and the result of each check
  00_setup.R              packages, balnet 0.0.4 from pkg/, sessionInfo.txt
  run_checks.R            the setup, tests, summaries, render and manifest in turn, output in logs/
  verify_summaries.R      recompute the Study 1 section summaries and the coverage table
                          from the raw results
  verify_render.R         render the thesis from the stored summaries in a clean copy
  smoke_test.R            rerun 15 replications of the first and last cell of every grid, n up
                          to 50 000, and compare them with the stored results
  smoke_report.csv        the smoke run on macOS, one row per stored file
  make_manifest.R         every file the thesis depends on, with sizes and md5 sums
  verification/
    displays.csv          every table and figure: page, qmd chunk, summaries, results files, runs
    compare_render.py     the render against the stored PDF, number by number and word by word
    make_verification.R   VERIFICATION.md from the logs
  R/
    dgp.R                 the design (Methods 3.2)
    estimators.R          the Study 1 estimators (Methods 3.3)
    selectors.R           the Study 2 selectors and lambda grids (Methods 3.4, 3.5.1)
    measures.R            bias, EmpSE, RMSE, MCSE, coverage at the empirical SE, R-bar, s-bar,
                          the two Study 2 rules (Methods 3.5)
    simulate.R            seed streams and the loop over a grid
    wong_chan.R           the Wong and Chan check (Appendix H)
  runs/
    study1.R              every Study 1 grid, one run each (Table 2)
    study2.R              every Study 2 grid, one run each (Table 3)
    wong_chan.R           the Wong and Chan cells
  tests/
    run_tests.R           runs the nine test files
    make_fixtures.R       cuts the rows the tests use from results/ into fixtures/
    fixtures/             those rows, so the tests run without results/
```

`replication/R/` rewrites the code that produced the results. That code stays in
`R/`, `runs/` and `handoff/` as the record of what ran, and the tests check that the
two draw the same data and return the same numbers. Results files name the
estimators `balnet0`, `glmnetcv_hajek`, `abw_balnet0`, `aipw_glmnetcv` and
`oracle_hajek`. `estimator_labels` in `replication/R/estimators.R` gives the printed
names.

## Three ways to reproduce

`main.R` runs these in order and prints the replications each grid used.

The thesis (`thesis/AAAAA Working Thesis.qmd`) renders from stored summaries in
`output/thesis_summaries/`, built from the raw results by
`handoff/scripts/build_summaries.R`. The coverage table (`coverage_empse`) is built
once by the qmd chunk `cov-data` instead. With the stored files in place no chunk
reads a raw results file.

1. Render the thesis from the stored summaries. Needs R, Quarto, a TeX installation with
   XeLaTeX, the Times New Roman font and the packages in `replication/00_setup.R`, not
   the raw results. `Rscript replication/verify_render.R` renders it in a temporary copy
   and leaves the stored PDF alone. `quarto render` from the project root writes to
   `output/rendered/` and replaces the stored PDF.
2. Recompute the summaries from the raw results. `Rscript replication/verify_summaries.R`
   prints the largest difference between recomputed and stored numbers for every
   Study 1 section summary and the coverage table. `handoff/scripts/build_summaries.R`
   rebuilds them all. It also checks the appendix entries, the Study 2 summaries and
   `wc_compare` by a second computation (`output/thesis_summaries/verification.md`).
3. Rerun the simulations. Weeks of machine time in total, see the run scripts.
   Every run uses L'Ecuyer-CMRG streams from a fixed base seed, so a rerun gives
   the stored draws exactly. Study 2 also fixes each cell's lambda grid from a pilot
   draw on Mersenne-Twister seed 1000 + i. The base seeds are in
   `replication/RESULTS_INDEX.md`.
   Each grid writes one file to `results/study1/` or `results/study2/`, and
   `replication/RESULTS_INDEX.md` lists the stored files each one reproduces. A smoke
   run of a few replications (`N_SIM=2 OUT_DIR=<folder>`) pairs with the stored files
   replication by replication. The `VAR=value Rscript` forms need a POSIX shell, such
   as Git Bash on Windows.

## Requirements

R 4.5 or 4.6. balnet 0.0.4 (the tarball in `pkg/` is the version that ran and needs
a C++17 compiler), Rcpp, RcppEigen, glmnet 5.0, data.table with R.utils, dplyr,
tidyr, purrr, tibble, here, testthat, withr, future and furrr. The render also needs
ggplot2, gt, ggh4x, knitr, rmarkdown, Quarto, a TeX installation with XeLaTeX and
the Times New Roman font, and the rebuild of the summaries needs readr.
`replication/00_setup.R` checks and installs the R packages. The stored results were
made under R 4.5.2 on Windows 11 and R 4.6.1 on macOS, and the Wong and Chan cells
under glmnet 4.1.10.

The last checks ran with these versions. The full list is in
`replication/sessionInfo.txt`.

| Software | Version |
|---|---------|
| R | 4.6.1 on macOS 26.6.2 |
| R packages | balnet 0.0.4, data.table 1.18.4, dplyr 1.2.1, furrr 0.4.0, future 1.75.0, ggh4x 0.3.1, ggplot2 4.0.3, glmnet 5.1, gt 1.3.0, here 1.0.2, knitr 1.51, Matrix 1.7-5, purrr 1.2.2, R.utils 2.13.0, readr 2.2.0, rmarkdown 2.31, testthat 3.3.2, tibble 3.3.1, tidyr 1.3.2, withr 3.0.3 |
| Quarto | 1.9.38, from RStudio |
| TeX | XeLaTeX, TeX Live 2022 |
| PDF comparison | Python 3.10, pdftotext 22.02 |

## Checks

```
Rscript main.R                             # setup tests summaries render and reruns. options at its top
Rscript replication/run_checks.R           # the next five in turn, needs results/
Rscript replication/00_setup.R
Rscript replication/tests/run_tests.R      # about 5 minutes
Rscript replication/verify_summaries.R     # about 7 minutes, needs results/
Rscript replication/verify_render.R        # needs Quarto and XeLaTeX
Rscript replication/make_manifest.R
python3 replication/verification/compare_render.py     # needs pdftotext (poppler)
Rscript replication/verification/make_verification.R   # writes VERIFICATION.md
Rscript replication/smoke_test.R           # about 20 minutes on 17 workers, needs results/
```

The times are from the last runs, on a MacBook Pro with an Apple M5 Pro chip (18
cores, 24 GB of memory) under macOS 26.6.2 (`replication/logs/machine.txt`).
`replication/logs/run_checks.txt` has the time of each step.

Without `results/` the tests read the rows they use from
`replication/tests/fixtures/`, and only the three recomputations in
`test_measures.R` skip.

The n = 1 million cell is not rerun by `smoke_test.R`, and
`replication/tests/test_runs.R` regenerates its first draw. Where a `balnet` path
stops short of exact balance, the fits just before its end are sensitive to
rounding, so a rerun on another platform can differ there, with the path ending at
the same lambda.

A rebuild of the summaries on macOS matches every number the displays read to 2e-9.
Two things differ. `app_p_axis` was built before the path endpoints `lam_end1` and
`lam_end0` were added to the p/n grid's file list, so a rebuild adds them and
redraws the MAE bootstrap, which no display prints. In Study 2 the grid point at
lambda = 1e-4 can be the same fit as the path endpoint, and on macOS the second
computation in `handoff/R/verify.R` breaks that tie differently in six cells. Three
summaries then fail one check and `wc_compare` is not stored.

## Not in this package

The raw results (`results/`, 7 GB). The original code of designs the thesis does not
use (an earlier IPW study, an earlier cross-validation design and exploratory Wong
and Chan runs). The treated proportion and path endpoint grids were run but are not
reported, and are not in the run scripts. Their files and runners remain in
`results/aipw` and `runs/aipw`.

## Licence

The code is under the BSD 3-Clause licence and the stored summaries and raw results
under CC BY 4.0. The thesis text is the author's, all rights reserved. balnet in
`pkg/` keeps its own MIT licence. The terms are in `LICENSE.txt`.

## References

Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths for
generalized linear models via coordinate descent. Journal of Statistical Software
33(1), 1-22. https://doi.org/10.18637/jss.v033.i01

R Core Team (2025). R: a language and environment for statistical computing.
R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Sverdrup, E. and Hastie, T. (2026). balnet: pathwise estimation of covariate
balancing propensity scores. https://doi.org/10.48550/arXiv.2602.18577

Wong, R. K. W. and Chan, K. C. G. (2018). Kernel-based covariate functional
balancing for observational studies. Biometrika 105(1), 199-213.
https://doi.org/10.1093/biomet/asx069
