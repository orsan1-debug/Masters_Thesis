# every study 1 grid (Table 2) to results/study1/<grid>.csv.gz. base seed 202
#   Rscript replication/runs/study1.R   every grid not yet run. about 5000 core-hours in all

source(here::here("replication", "R", "dgp.R"))
source(here::here("replication", "R", "estimators.R"))
source(here::here("replication", "R", "simulate.R"))

out_dir <- Sys.getenv("OUT_DIR", here::here("results", "study1"))
workers <- as.integer(Sys.getenv("WORKERS", max(1L, parallel::detectCores() - 1L)))
n_sim   <- as.integer(Sys.getenv("N_SIM", NA))

# 4 confounders in p = 50 w/ Kang-Schafer misspec. denser design w/ Tan's transform
four_confounders <- function(cell)
  generate_sample(n = cell$n, p = if (is.null(cell$p)) 50 else cell$p, s = 4, c = cell$overlap,
                  kappa = if (is.null(cell$strength)) 1 else cell$strength,
                  misspecification = if (cell$spec == "misspecified") "kang_schafer" else "none")
denser <- function(cell)
  generate_sample(n = cell$n, p = cell$p, s = cell$s, nu = cell$decay,
                  covariance = c(iid = "independent", ar1 = "ar1", equi = "equicorrelated")[[cell$covcor]],
                  outcomes = c("linear", "quadratic"), loadings = "all_active",
                  misspecification = if (cell$spec == "misspecified") "tan" else "none")

both_specs <- c("correct", "misspecified")   # share a stream block so same draws

denser_points <- function(covariance) {
  g <- expand.grid(covcor = covariance, decay = c(1, 0.5), s = c(4, 16, 64), pn = c(0.1, 0.5, 1),
                   n = c(500, 1000), stringsAsFactors = FALSE)
  g$p <- as.integer(g$pn * g$n)
  g <- g[g$s <= g$p, c("n", "pn", "p", "s", "decay", "covcor")]
  g$block <- seq_len(nrow(g))
  g
}
levels_c <- c(0.05, 0.15, 0.25, 0.5, 0.75, 1, 2)
n_path   <- c(250, 500, 1000, seq(2000, 20000, by = 1000), 30000, 40000, 50000, 75000, 1e5, 1.5e5, 2e5)

# blocks follow the stored runs so draws and folds match. stored = results/aipw files per spec
grids <- list(
  # misspec stops at n = 50k like the stored runs
  sample_size = list(draw = four_confounders, replications = 10000,
    cells = subset(under_each(data.frame(n = c(500, 1000, 5000, 10000, 50000, 1e5, 2e5, 4e5), overlap = 1),
                              spec = both_specs),
                   spec == "correct" | n <= 50000),
    stored = list(correct = c("aipw_correct_n400k_dgp2", "aipw_correct_n_dgp2", "aipw_add_correct_n400k_dgp2"),
                  misspecified = c("aipw_misspec_n_dgp2", "aipw_add_misspec_n_dgp2"))),

  # n = 1 million was its own run (10 batches of 1000) so it starts at block 1
  sample_size_path = list(draw = four_confounders, replications = 10000,
    cells = under_each(data.frame(n = c(n_path, 1e6), overlap = 1, block = c(seq_along(n_path), 1)),
                       spec = "misspecified"),
    stored = list(misspecified = c("aipw_full_misspec_nfine_dgp2", "aipw_misspec_n_dgp2n1mill"))),

  # n = 1000 and 10k were separate runs so each starts at block 1
  overlap = list(draw = four_confounders, replications = 10000,
    cells = under_each(data.frame(n = rep(c(1000, 10000), each = 7), overlap = levels_c, block = rep(1:7, 2)),
                       spec = both_specs),
    stored = list(correct = c("aipw_correct_overlap7_dgp2", "aipw_diag_correct_overlap7_dgp2",
                              "aipw_correct_overlap8_dgp2"),
                  misspecified = c("aipw_misspec_overlap7_dgp2", "aipw_full_misspec_overlap7_dgp2",
                                   "aipw_misspec_overlap8_dgp2"))),

  n_by_overlap = list(draw = four_confounders, replications = 1000,
    cells = under_each(expand.grid(n = round(signif(exp(seq(log(100), log(5000), length.out = 12)), 2)),
                                   overlap = signif(10^seq(-2, 0, by = 1 / 6), 3)), spec = both_specs),
    stored = list(correct = "aipw_correct_n_by_overlap_dgp2", misspecified = "aipw_misspec_n_by_overlap_dgp2")),

  dimension = list(draw = four_confounders, replications = 10000,
    cells = under_each(local({ g <- expand.grid(n = c(500, 1000), pn = c(0.05, 0.1, 0.5, 1, 2, 4))
                               data.frame(n = g$n, p = as.integer(g$pn * g$n), overlap = 1) }),
                       spec = both_specs),
    stored = list(correct = c("aipw_correct_p_dgp2", "aipw_diag_correct_p_dgp2"),
                  misspecified = c("aipw_misspec_p_dgp2", "aipw_diag_misspec_p_dgp2"))),

  outcome_signal = list(draw = four_confounders, replications = 10000,
    cells = under_each(data.frame(n = 1000, overlap = 1, strength = c(0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1)),
                       spec = both_specs),
    stored = list(correct = "aipw_correct_strength_dgp2", misspecified = "aipw_misspec_strength_dgp2")),

  # iid + ar1 were one run and equi another. each starts at block 1
  denser = list(draw = denser, replications = 1000,
    cells = under_each(rbind(denser_points(c("iid", "ar1")), denser_points("equi")), spec = both_specs),
    stored = list(correct = c("aipw_correct_sparsity_hd_dgp2", "aipw_correct_sparsity_hd_equi_dgp3"),
                  misspecified = c("aipw_misspec_tanobs_sparsity_hd_dgp3",
                                   "aipw_misspec_tanobs_sparsity_hd_equi_dgp3"))))

#' one grid into results/study1. skips if its file exists
run_study1_grid <- function(name) {
  g <- grids[[name]]
  out <- file.path(out_dir, paste0(name, ".csv.gz"))
  if (file.exists(out)) return(message("skip ", basename(out)))
  run_grid(g$draw, g$cells, replications = if (is.na(n_sim)) g$replications else min(n_sim, g$replications),
           base_seed = 202, out_file = out, estimate = estimate_study1_with_reference, workers = workers,
           stream_block = g$replications, cells = cells_to_run(g$cells))
}

# only runs under Rscript. main.R smoke_test.R and the tests source it for grids
if (sys.nframe() == 0L)
  for (name in grep(Sys.getenv("GRIDS", "."), names(grids), value = TRUE)) run_study1_grid(name)
