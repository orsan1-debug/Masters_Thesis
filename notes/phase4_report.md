# Phase 4 verification report

Repository state verified: commit `1040a58` (Phase 3b). Date: 2026-09-11. Nothing was fixed; findings only.
Verification byproducts (Quarto's per-directory `.gitignore` stubs) were removed; see section 5.

## 1. Sourcing R/*.R in a fresh session

Order: `dgp.R`, `estimators_cv.R`, `estimators_ipw.R`, `registry.R`, `simulate.R`, `summarise.R`, `plots.R`
(alphabetical, with `summarise.R` before `plots.R`). Run with `Rscript --vanilla` from the repo root.

- Errors: none. All seven files source; 57 functions and 17 non-function objects are defined afterwards.
- Warnings: "package 'X' was built under R version 4.5.3" for sbw, slam, dplyr, tidyr, ggplot2, patchwork, gt
  (the installed R is 4.5.2). Masking messages from dplyr (filter, lag, intersect, setdiff, setequal, union;
  `summarize` from sbw) and tidyr (expand, pack, unpack from Matrix).
- Side effects (top-level, not inside functions):
  - `estimators_cv.R`: `library(sbw)` attaches sbw, Matrix, quadprog, slam; `source(here::here("R", "dgp.R"))`
    re-sources dgp.R (harmless, redefines the same functions).
  - `summarise.R`: attaches dplyr, tidyr, ggplot2, patchwork, gt; assigns `select`/`filter` aliases into the
    global environment.
  - options changed by package loads: `ambiguousMethodSelection` (Matrix), `gt.*` (3) and
    `callr.condition_handler_cli_message` (gt/cli). None set by project code.
  - Files created or modified: none. Working directory: unchanged. `library()` inside a function: none.
- Observation: no R/ file attaches balnet or glmnet; `estimators_ipw.R` calls `cv.balnet()`, `balweights()`
  and `cv.glmnet()` unqualified, so a caller must attach them (the ipw run scripts rely on `simulate_grid()`
  passing `packages = c("balnet", "glmnet")` to the workers, and on the session having them attached).

## 2. runs/ smoke runs (2 reps)

Not run. In every runs/ script the replication count is a hard-coded literal, so per instruction the line is
listed and the run stopped. No file in results/ or registry.csv was touched.

| subdir | smallest script | hard-coded count |
|---|---|---|
| runs/ipw | `correctspecDGP1.R` (827 bytes) | line 16: `num_sim   = 1000,` inside the `run_batch()` call |
| runs/cv/basic_dgp | `run_bad_overlap_noise10_1000reps__tune4_19.R` (2,491 bytes) | line 20: `n_rep <- 1000` |
| runs/cv/wz_replication | `run_wc_path.R` (4,767 bytes) | line 23: `n_rep       <- 1000` |

All 13 ipw scripts pass `num_sim = 1000` literally to `run_batch()`; all 19 basic_dgp and 5 wz scripts set
`n_rep <-` literally (50 to 1000). Registry behaviour, from the code: every ipw script goes through
`run_batch()`, whose `on.exit()` calls `register_run()`, which appends a row to `registry.csv`
(default `here::here("registry.csv")`) on every run, including failed ones. The cv scripts never touch
registry.csv. Output paths: ipw scripts write `here::here("results", "ipw", <batch>.csv.gz)` plus the
`_cells/` checkpoint directory and `_cvcurves.rds`/`_session.txt` sidecars; cv scripts write to
`res_dir = here::here("results", "cv", ...)`. A 2-rep smoke run therefore needs the count and `out_file`
(ipw) or `res_dir` (cv) made overridable before it can be pointed at a temp directory.

## 3. Quarto renders (one per analysis/ subdirectory, output to temp directories)

## 3a. Render: analysis/ipw

- **File rendered:** `analysis/ipw/SNR.qmd` (59 lines; smallest of 11 .qmd in the dir). Sources `R/summarise.R`, `R/plots.R` and reads `results/ipw/snrDGP2correct.csv.gz` (4.9 MB) via `here::here()`.
- **Command:** `cd "<repo root>" && quarto render "analysis/ipw/SNR.qmd" --output-dir "<scratch>/out" 2>&1 | tee "<scratch>/render.log"`
- **Runtime:** start 17:39:49, end 17:39:59 (AUSEST) -> ~10 s wall clock. Exit code 0.
- **Result:** OK. Both chunks (`unnamed-chunk-1` setup, `unnamed-chunk-2` plots/table) executed; knitr -> pandoc -> HTML with no failing chunk.
- **Output produced:** `<scratch>/out/SNR.html` (53,784 bytes) plus `<scratch>/out/SNR_files/` containing 7 figure PNGs (`figure-html/unnamed-chunk-2-{1..7}.png`, 10-28 KB each) and the standard quarto/bootstrap libs. HTML references all 7 images and contains 1 table (the `perf_table`). Note: `here::here()` resolution worked from the repo root, and nothing in `results/` was written.
- **Repo byproducts found (git status --short, before vs after):**
  - `?? analysis/ipw/.gitignore` (30 bytes, contents `/.quarto/` + `**/*.quarto_ipynb`; mtime 17:39:50, i.e. inside this render's window, in the rendered file's directory; absent from the before-status and a full before-render file-tree snapshot). Proven byproduct of this render: quarto writes this when rendering a standalone .qmd inside a git repo with no `_quarto.yml`. **Attempted `rm` but the action was blocked by the auto-mode permission classifier, so it was NOT deleted and remains in the repo.** It is untracked, harmless, and safe to delete.
  - `?? analysis/cv/basic_dgp/.gitignore` (mtime 17:39:46, i.e. 3 s BEFORE this render started) and `?? analysis/cv/wz_replication/.gitignore` (mtime 17:39:52) - identical 30-byte contents. Both are in directories this render never touched; they appeared during the same seconds, so they almost certainly come from concurrent renders of those directories (sibling verification tasks), not from this one. Not touched.
  - No `.quarto/` directory, no `SNR.html`, no `SNR_files/` were created next to the qmd or at the repo root (verified by diffing a `find` listing of the whole tree before/after; the only tree diff is the three `.gitignore` files above). No pre-existing tracked file was modified.
- **R warnings during render:** none printed. The setup chunk sets `message: false` / `warning: false` so package/`load_sim` messages would be suppressed by design; the second (unsuppressed) chunk emitted no warnings, `render.log` contains no "warning"/"error" lines, and the HTML contains 0 occurrences of "Warning".
- **Observations (not errors):** the YAML title of `SNR.qmd` is "Overlap (DGP2)" although the file is about SNR - likely a copy-paste leftover. No `_quarto.yml` exists anywhere in the repo, so every standalone `quarto render` will keep (re)creating a per-directory `.gitignore`; adding `.quarto/` to the root `.gitignore` or a root `_quarto.yml` would stop that.

## 3b. Render: analysis/cv/basic_dgp

- **File rendered:** `analysis/cv/basic_dgp/cv_summary.qmd` (default format: pdf, pdf-engine xelatex, mainfont Times New Roman)
- **Command:** `cd "<repo root>" && quarto render "analysis/cv/basic_dgp/cv_summary.qmd" --output-dir "<scratch dir>/out" 2>&1 | tee "<scratch dir>/render.log"`
- **Runtime:** ~9 s wall clock (start 17:39:45, end 17:39:54 AUSEST)
- **Result:** OK
- **Failing chunk / error:** none. All 9 knitr steps completed (`setup`, `unnamed-chunk-1..3`); pandoc -> latex; xelatex ran twice (XeTeX 3.141592653-2.6-0.999998, TeX Live 2026). TeX engine present; all csv inputs (via `here::here("output","cv","basic_dgp","summaries")`) and png figures (`../../../output/cv/basic_dgp/figures/`) resolved.
- **Output produced:** `<scratch dir>/out/cv_summary.pdf`, 499,674 bytes (488 KB). Log at `<scratch dir>/render.log`.
- **R warnings/messages printed during render:** none (log contains no "warning" or "error" lines; document sets `warning: false`, `message: false`, so R-level warnings would be suppressed anyway).
- **Repo byproducts (git status --short before vs after):**
  - `?? analysis/cv/basic_dgp/.gitignore` (30 bytes, `/.quarto/` + `**/*.quarto_ipynb`, created 17:39:46 in the render directory) -- created by this render. **Deleted** it; final `git status --short` is identical to the pre-render snapshot for this path.
  - `?? analysis/cv/wz_replication/.gitignore` (created 17:39:52) and `?? analysis/ipw/.gitignore` (created 17:39:50), same 30-byte Quarto content. These appeared during the render window but sit in directories this render never touched; consistent with concurrent renders by parallel verification tasks. **Not deleted** (cannot prove this render created them). By the final status check `wz_replication/.gitignore` had already disappeared on its own (removed by whichever task owns it), confirming the concurrent-render explanation; only `analysis/ipw/.gitignore` remained at last check and should be handled by the ipw render task.
  - No `.quarto/`, `*_files/`, `.pdf`, `.tex`, `.log`, or `.knit.md` left next to the qmd (intermediates were cleaned by quarto; the pdf went to the scratch out dir).
- **Pre-existing note (not a byproduct):** `output/cv/basic_dgp/cv_summary.pdf` already existed before the render and is git-ignored (`output/` is `!!`); it was not touched.
- **Pre-render untracked entries (unchanged, for reference):** 15 `??` entries, all under `archive/` (e.g. `archive/exploring_cv/cv_summary.{aux,log,tex}`, `archive/extension_cv/`).

## 3c. Render: analysis/cv/wz_replication

- **File rendered:** `analysis/cv/wz_replication/wc_findings.qmd` (the only .qmd in that directory; YAML: `format: pdf`, `pdf-engine: xelatex`, `mainfont: Times New Roman`; rendered with its default format, unchanged).
- **Command:** `cd "<repo root>" && quarto render "analysis/cv/wz_replication/wc_findings.qmd" --output-dir "<scratch>/out" 2>&1 | tee "<scratch>/render.log"`
- **Runtime:** ~7 s wall clock (start 17:39:52, end 17:39:59 AUSEST, 2026-09-11). Quarto exit code 0.
- **Result:** **OK**
- **Failing chunk / error:** none. All 7 knitr units ran (labelled chunks: `setup`, `sbw`, `att`); knit -> pandoc -> xelatex (2 passes, TeX Live 2026) completed. No missing TeX engine, no missing data file (`here::here()` resolved `wc_grids.R`, `wc_sbw_v2.csv.gz`, `wc_att_v2.rds`, and the `wc_noise_*`, `wc_on_*`, `wc_basis_*` .rds grids under `results/cv/wz_replication/`).
- **Output produced:** `<scratch>/out/wc_findings.pdf`, 77,414 bytes (~76 KB). Full log at `<scratch>/render.log`.
- **Repo byproducts:** `git status --short` before = 15 entries, after = 18. The 3 new entries were all untracked `.gitignore` files (content `/.quarto/` + `**/*.quarto_ipynb`, the stub Quarto writes next to a rendered doc):
  - `analysis/cv/wz_replication/.gitignore` (mtime 17:39:52.69, inside this render window, in the rendered doc's directory) -- **created by this render; deleted.** Directory now shows no status entries.
  - `analysis/cv/basic_dgp/.gitignore` (mtime 17:39:46) and `analysis/ipw/.gitignore` (mtime 17:39:50) -- created *before* this render started, evidently by concurrent renders of other directories. **Not touched**; flagged for whoever owns those renders.
  - No `.quarto/`, `*_files/`, `.pdf`, `.tex`, `.log`, or `.knit.md` left in the repo (intermediates went to the scratch output dir / were cleaned by Quarto). No pre-existing entries changed.
- **R warnings during render:** none printed (the document sets `execute: warning: false`, so chunk-level R warnings would be suppressed; nothing appeared from knitr, pandoc, or xelatex either).

## 4. Static checks (non-archived .R and .qmd under R/, runs/, analysis/, tests/, thesis/)

- Absolute paths or `setwd()`: 0 (comments and URLs included in the scan).
- `library()` / `require()` inside a function body: 0.
- Roxygen coverage in R/: 54 top-level functions, 49 without a `#'` block. The 5 that have one are
  `fit_sbw`, `cstat`, `estimate_wc_sbw` (estimators_cv.R) and `gen_data`, `run_par` (dgp.R). The 49 are
  listed in the Phase 4 log; fixing is Phase 5.
- Parse: every non-archived file parses. The two archived files that do not (`archive/tests/balnet_fix.R`,
  `archive/analysis/cv_tuning.qmd`) were broken before the restructure.

## 5. Git status

`git status --short` after verification: 15 untracked files, all under `archive/` and all untracked before
the restructure (`archive/crit_paths_seed.png`, `archive/cv_extension/files_1.zip`,
`archive/exploring_cv/{balnet_scoring_functions.R, balnet_selection_functions.R, cv_flip.csv, cv_summary.aux,
cv_summary.log, cv_summary.tex, cv_summary_working.qmd, factors_paths.csv, images/clipboard-*.png x4}`,
`archive/extension_cv/`). No tracked file was modified by the verification.

Each `quarto render` created a 30-byte `.gitignore` stub (`/.quarto/`, `**/*.quarto_ipynb`) next to the
rendered document; all three stubs were removed after the renders (two by the render tasks, one by the
coordinator). Because the repo has no `_quarto.yml`, every future standalone render will recreate them;
adding `.quarto/` to the root `.gitignore` or a root `_quarto.yml` would stop that (Phase 5 candidate).
Also noted by the ipw render: `analysis/ipw/SNR.qmd` carries the YAML title "Overlap (DGP2)".
