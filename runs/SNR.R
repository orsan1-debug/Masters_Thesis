## SNR GRID ---------------------------------------------------------------
## DGP2, correct spec, iid, p = 50, overlap = 1, vary signal strength.
## SNR = strength^2 * Var(f); at strength = 1, linear = 1.75, quad1 = 5.
source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")
grid <- expand.grid(n = c(1000, 10000), overlap = 1,
                    strength = c(0.25, 0.5, 1, 2, 4))
dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap,
       strength = cell$strength)
res_snr <- run_batch(dgp_gen(FALSE), grid,
                     num_sim   = 1000,
                     base_seed = 203,
                     out_file  = "results/snrDGP2correct.csv.gz",
                     meta = list(label = "DGP2snrcorrectspec"))