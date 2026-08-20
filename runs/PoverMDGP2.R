POVERM DGP2 ----------------------------------------------------------------
  ## DGP2, correct spec, iid, p >> n grid, overlap = 1, 1000 reps.
  ## Extends seed-204 dimensionality ratio axis {0.1-2} to {4, 8} at same n.
  source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")
grid   <- expand.grid(n = c(500, 1000),
                      pn_ratio = c(4, 8))
grid$p <- grid$n * grid$pn_ratio
dgp_gen <- function(cell)
  dgp2(n = cell$n, p = cell$p, s = 4, signs = "pos",
       outcome = c("linear", "quad1"),
       covcor = "iid", misspec = FALSE, overlap = 1)
res_pm <- run_batch(dgp_gen, grid,
                    num_sim   = 1000,
                    base_seed = 205,
                    out_file  = "results/PoverMDGP2.csv.gz",
                    meta = list(label = "DGP2PoverM"))

