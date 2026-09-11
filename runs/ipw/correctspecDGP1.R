## CORRECT SPEC DGP1 ---------------------------------------------------------
## DGP1, correct spec, iid, p = 50, overlap = 1, 1000 reps.
## Supersedes legacy E1 + REDUX. Floor 1e-4 via estimate_all defaults.
source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))
n_sim   <- as.integer(Sys.getenv("N_SIM", "1000"))          # N_SIM=2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "ipw"))  # OUT_DIR=<tmp> keeps results/ untouched

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000), overlap = 1)

dgp_gen_e1 <- function(cell)
  dgp1(n = cell$n, p = 50, outcome = c("linear", "quad1", "exp"),
       overlap = cell$overlap)

res <- run_batch(dgp_gen_e1, grid,
                 num_sim   = n_sim,
                 base_seed = 101,
                 out_file  = file.path(out_dir, "correctspecDGP1.csv.gz"),
                 dgp = "dgp1", script = "runs/ipw/correctspecDGP1.R")
