# every study 2 grid (Table 3) to results/study2/<grid>.csv.gz + <grid>_lamgrid.csv
#   Rscript replication/runs/study2.R   overlap ~2 h on 8 cores. penalty_by_overlap most of a day

source(here::here("replication", "R", "dgp.R"))
source(here::here("replication", "R", "estimators.R"))
source(here::here("replication", "R", "selectors.R"))
source(here::here("replication", "R", "simulate.R"))

out_dir <- Sys.getenv("OUT_DIR", here::here("results", "study2"))
workers <- as.integer(Sys.getenv("WORKERS", max(1L, parallel::detectCores() - 1L)))
n_sim   <- as.integer(Sys.getenv("N_SIM", NA))

# own_strength = one strength per cell. shared = all 7 strengths on one draw via share_strengths()
# dense = the denser design w/ Tan misspec. others default to n = 1000 p = 50 and iid covs
own_strength <- function(cell)
  generate_sample(n = 1000, p = 50, s = 4, c = cell$overlap, kappa = cell$strength,
                  misspecification = if (identical(cell$spec, "misspecified")) "kang_schafer" else "none")
shared <- function(cell)
  share_strengths(generate_sample(n = if (is.null(cell$n)) 1000 else cell$n,
                                  p = if (is.null(cell$p)) 50 else cell$p, s = 4, c = cell$overlap,
                                  misspecification = if (identical(cell$spec, "misspecified")) "kang_schafer" else "none"))
dense <- function(cell)
  share_strengths(generate_sample(n = cell$n, p = cell$p, s = cell$s, nu = cell$decay, c = cell$overlap,
                                  covariance = c(iid = "independent", ar1 = "ar1", equi = "equicorrelated")[[cell$covcor]],
                                  loadings = "all_active",
                                  misspecification = if (identical(cell$spec, "misspecified")) "tan" else "none"))

strengths  <- c(1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0)   # the noise grids
overlaps   <- c(1, 0.75, 0.5, 0.25, 0.15, 0.05)
with_floor <- function(cells) { cells$min_ratio <- ifelse(cells$alpha == 0, 1e-4, 1e-2); cells }
dense_points <- expand.grid(s = c(4, 16, 64), covcor = c("iid", "ar1", "equi"), stringsAsFactors = FALSE)
dense_points$n <- 1000; dense_points$p <- 100; dense_points$decay <- 0.5; dense_points$overlap <- 0.25

# 2500 reps unless given. each grid keeps the seed and reps of its stored run
# stored = files in results/cv/our dgp/r2500 (r1000 at 1000 reps)
grids <- list(
  overlap = list(draw = own_strength, seed = 301,
    cells = data.frame(overlap = overlaps, strength = 1),
    stored = "cv_overlap_dgp2"),
  noise = list(draw = own_strength, seed = 302,
    cells = data.frame(overlap = 1, strength = strengths),
    stored = "cv_noise_dgp2"),
  noise_by_overlap = list(draw = own_strength, seed = 303,
    cells = expand.grid(overlap = c(1, 0.5, 0.25, 0.15), strength = strengths),
    stored = "cv_noise_x_overlap_dgp2"),
  penalty_by_overlap = list(draw = shared, seed = 305,
    cells = with_floor(expand.grid(alpha = c(0, 0.25, 0.5, 0.75, 1), overlap = c(1, 0.5, 0.25, 0.15))),
    stored = "cv_alpha_x_overlap_dgp2"),
  sample_size = list(draw = shared, seed = 306,
    cells = with_floor(expand.grid(n = c(250, 500, 1000, 2000), overlap = c(0.25, 0.15), alpha = c(0.5, 1))),
    stored = "cv_n_dgp2"),
  dimension = list(draw = shared, seed = 307,
    cells = with_floor(expand.grid(p = c(50, 250, 500, 1000), overlap = c(0.5, 0.25), alpha = c(0.5, 1))),
    stored = "cv_p_dgp2"),
  misspecification = list(draw = shared, seed = 308,
    cells = data.frame(overlap = overlaps, spec = "misspecified"),
    stored = "cv_misspec_dgp2"),
  # not reported. the 3 fold counts were separate runs on one seed so they share a block
  folds = list(draw = shared, seed = 309,
    cells = under_each(data.frame(overlap = c(1, 0.5, 0.25, 0.15)), folds = c(10, 5, 3)),
    stored = list(`10` = "cv_folds_k10_dgp2", `5` = "cv_folds_k5_dgp2", `3` = "cv_folds_k3_dgp2")),
  # n = 5000 and 10k were their own runs w/ their own seeds
  large_n_5000 = list(draw = shared, seed = 310,
    cells = data.frame(n = 5000, overlap = c(0.25, 0.15)),
    stored = "cv_n5000_dgp2"),
  large_n_10000 = list(draw = shared, seed = 311, replications = 1000,
    cells = data.frame(n = 10000, overlap = c(0.25, 0.15)),
    stored = "cv_n10000_dgp2"),
  # the 2 specs were separate runs on one seed so they share a block
  denser = list(draw = dense, seed = 312,
    cells = under_each(dense_points, spec = c("correct", "misspecified")),
    stored = list(correct = "cv_sparsity_dgp3", misspecified = "cv_sparsity_misspec_dgp3")))

#' one grid into results/study2. skips if its file exists
run_study2 <- function(name) {
  g <- grids[[name]]
  replications <- if (is.null(g$replications)) 2500 else g$replications
  run_study2_grid(g$draw, g$cells, replications = if (is.na(n_sim)) replications else min(n_sim, replications),
                  base_seed = g$seed, out_file = file.path(out_dir, paste0(name, ".csv.gz")),
                  workers = workers, stream_block = replications, cells = cells_to_run(g$cells))
}

# only runs under Rscript. smoke_test.R and tests/test_runs.R source it for grids
if (sys.nframe() == 0L)
  for (name in grep(Sys.getenv("GRIDS", "."), names(grids), value = TRUE)) run_study2(name)
