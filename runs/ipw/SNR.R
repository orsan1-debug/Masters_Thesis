## SNR GRID ---------------------------------------------------------------
## DGP2, correct spec, iid, p = 50, overlap = 1, vary signal strength.
## SNR = strength^2 * Var(f); at strength = 1, linear = 1.75, quad1 = 5.
source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))
n_sim   <- as.integer(Sys.getenv("N_SIM", "1000"))          # N_SIM=2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "ipw"))  # OUT_DIR=<tmp> keeps results/ untouched
grid <- expand.grid(n = c(1000, 10000), overlap = 1,
                    strength = c(0.25, 0.5, 1, 2, 4))
dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap,
       strength = cell$strength)
res_snr <- run_batch(dgp_gen(FALSE), grid,
                     num_sim   = n_sim,
                     base_seed = 203,
                     out_file  = file.path(out_dir, "snrDGP2correct.csv.gz"),
                     dgp = "dgp2", script = "runs/ipw/SNR.R")