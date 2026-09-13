## tests/check_aipw_results.R -------------------------------------------------
## Acceptance checks for one AIPW result file (Morris et al. 2019 s.5): no
## failed reps, the 35 row names, estimator rows identical to results/aipw_10k,
## diagnostic ranges, summarise_sim clean; then medians by overlap.
## Smoke: 5 reps into tempdir(), cell1_only = TRUE. Re-run: new file, FALSE.

rm(list = ls())
source(here::here("R", "packages.R"))
source(here::here("R", "summarise.R"))

# *** Settings ***
result_file <- file.path(tempdir(), "pilot.csv.gz")
ref_file    <- file.path(here::here("results", "aipw_10k"),
                         "aipw_correct_overlap_dgp2.csv.gz")  # or _misspec_
cell1_only  <- TRUE   # below num_sim = 10,000 only cell 1 shares the
#   aipw_10k streams (simulate.R line 77)
p           <- 50     # not stored in the file

# *** Expected rows ***
est  <- c("balnet0", "glmnetcv_ht", "glmnetcv_hajek", "aipw_glmnetcv",
          "abw_balnet0", "oracle_ht", "oracle_hajek")
arms <- c("1", "0", "_glm1", "_glm0", "_true1", "_true0")
nnz  <- c("nnz_bal1", "nnz_bal0", "nnz_ps", "nnz_or1", "nnz_or0")
tim  <- c("time_bal", "time_ps", "time_or")
diag <- c(paste0(c("smd", "ess", "wmax"), rep(arms, each = 3)),
          "lam_end1", "lam_end0", nnz, tim)

# *** Data ***
d     <- load_sim(result_file)
n_sim <- max(d$sim)
cmp   <- d |> filter(estimator %in% est) |>
  inner_join(load_sim(ref_file) |> filter(sim <= n_sim, estimator %in% est),
             by = c("overlap", "sim", "estimator", "outcome"),
             suffix = c("_new", "_ref"))
id    <- if (cell1_only) filter(cmp, overlap == 0.25) else cmp  # grid row 1
dw    <- diag_wide(d)
const <- d |>
  filter(estimator %in% setdiff(diag, c("nnz_or1", "nnz_or0"))) |>
  group_by(overlap, sim, estimator) |>
  summarise(v = n_distinct(tau_hat), .groups = "drop")

# *** Checks ***
stopifnot(
  "no failed reps"              = all(d$err == ""),
  "35 row names, no extras"     = setequal(unique(d$estimator), c(est, diag)),
  "is_diag = the 28 diags"      =
    setequal(unique(d$estimator[is_diag(d$estimator)]), diag),
  "every estimator row matched" =
    nrow(cmp) == n_distinct(d$overlap) * n_sim * length(est) *
    n_distinct(d$outcome),
  "estimator rows equal ref"    =
    max(abs(id$tau_hat_new - id$tau_hat_ref)) < 1e-10,
  "scalar diags constant"       = all(const$v == 1),
  "1 <= ess <= n"               = all(dw[paste0("ess", arms)] >= 1 &
                                        dw[paste0("ess", arms)] <= max(d$n)),
  "wmax >= 1"                   = all(dw[paste0("wmax", arms)] >= 1),
  "smd >= 0"                    = all(dw[paste0("smd", arms)] >= 0),
  "0 <= nnz <= p"               = all(dw[nnz] >= 0 & dw[nnz] <= p),
  "times > 0"                   = all(dw[tim] > 0),
  "summarise_sim: 7 rows, n_na 0" = {
    s <- summarise_sim(d)
    setequal(unique(s$estimator), est) && all(s$n_na == 0)
  }
)
cat("all checks passed:", basename(result_file), "\n")

# *** Table ***
dw |> filter(outcome == "linear") |>
  group_by(overlap) |>
  summarise(across(c(smd1, smd_glm1, smd_true1, ess1, ess_glm1, ess_true1,
                     wmax1, wmax_glm1, wmax_true1, nnz_bal1, nnz_ps, nnz_or1,
                     time_bal, time_ps, time_or), median)) |>
  print(width = Inf)