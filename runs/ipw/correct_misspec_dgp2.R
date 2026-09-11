## CORRECT and MISSPEC SPEC DGP2 ---------------------------------------------------------
## DGP2, correct spec + MISSPEC, iid, p = 50, overlap = 1, 1000 reps.


source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))
n_sim   <- as.integer(Sys.getenv("N_SIM", "1000"))          # N_SIM=2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "ipw"))  # OUT_DIR=<tmp> keeps results/ untouched

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000), overlap = 1)

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

res_cor <- run_batch(dgp_gen(FALSE), grid,
                     num_sim   = n_sim,
                     base_seed = 202,
                     out_file  = file.path(out_dir, "correctspecDGP2.csv.gz"),
                     dgp = "dgp2", script = "runs/ipw/correct_misspec_dgp2.R")

res_mis <- run_batch(dgp_gen(TRUE), grid,
                     num_sim   = n_sim,
                     base_seed = 202,
                     out_file  = file.path(out_dir, "misspecDGP2.csv.gz"),
                     dgp = "dgp2", script = "runs/ipw/correct_misspec_dgp2.R")