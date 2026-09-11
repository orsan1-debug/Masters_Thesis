## HIGH-P SPARSITY DGP2 --------------------------------------------------------
source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))
n_sim   <- as.integer(Sys.getenv("N_SIM", "1000"))          # N_SIM=2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "ipw"))  # OUT_DIR=<tmp> keeps results/ untouched
grid <- expand.grid(n = 1000,
                    p = c(50, 250, 500),
                    s = c(4, 16, 64, 256),
                    decay = c(0.5, 1))
grid <- grid[grid$s <= grid$p, ]   
dgp_gen <- function(cell)
  dgp2(n = cell$n, p = cell$p, s = cell$s, signs = "pos",
       outcome = c("linear"), outcome_set = "track_s",
       decay_ps = cell$decay, decay_out = cell$decay,
       covcor = "iid", misspec = FALSE, overlap = 1)

res_sphp <- run_batch(dgp_gen, grid,
                      num_sim   = n_sim,
                      base_seed = 208,
                      out_file  = file.path(out_dir, "sparsityE6DGP2.csv.gz"),
                      dgp = "dgp2", script = "runs/ipw/sparsity_2.R")

dgp_gen_ar1 <- function(cell)
  dgp2(n = cell$n, p = cell$p, s = cell$s, signs = "pos",
       outcome = "linear", outcome_set = "track_s",
       decay_ps = cell$decay, decay_out = cell$decay,
       covcor = "ar1", misspec = FALSE, overlap = 1)
res_sphp_ar1 <- run_batch(dgp_gen_ar1, grid,
                          num_sim   = n_sim,
                          base_seed = 209,
                          out_file  = file.path(out_dir, "sparsityE6ar1DGP2.csv.gz"),
                          dgp = "dgp2", script = "runs/ipw/sparsity_2.R")