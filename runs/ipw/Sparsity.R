## SPARSITY DGP2 ---------------------------------------------------------------
## DGP2, correct spec, iid, p = 50, s x decay grid, overlap = 1, 1000 reps.
## s = 50 = p is the exactly dense case; decay 0.5 flattens the coefficient
## profile; track_s, outcome shares weights. Supersedes failed dgp2spars batch.
source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))
n_sim   <- as.integer(Sys.getenv("N_SIM", "1000"))          # N_SIM=2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "ipw"))  # OUT_DIR=<tmp> keeps results/ untouched
grid <- expand.grid(n = c(500, 1000, 5000, 10000),
                    s = c(4, 8, 16, 32, 50),
                    decay = c(1, 0.5),
                    covcor = c("iid", "ar1"),
                    stringsAsFactors = FALSE)
dgp_gen <- function(cell)
  dgp2(n = cell$n, p = 50, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1"),
       outcome_set = "track_s", decay_ps = cell$decay, decay_out = cell$decay,
       covcor = cell$covcor, misspec = FALSE, overlap = 1)
res_sp <- run_batch(dgp_gen, grid,
                    num_sim   = n_sim,
                    base_seed = 206,
                    out_file  = file.path(out_dir, "sparsityDGP2.csv.gz"),
                    dgp = "dgp2", script = "runs/ipw/Sparsity.R")




