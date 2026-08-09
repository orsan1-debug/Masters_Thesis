## DIMENSIONALITY DGP2 --------------------------------------------------------
## DGP2, correct + misspec, iid, archived E5 p/n grid, overlap = 1, 1000 reps.
source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")
grid   <- expand.grid(n = c(500, 1000),
                      pn_ratio = c(0.1, 0.5, 1, 2))
grid$p <- grid$n * grid$pn_ratio
dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = cell$p, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = 1)
res_dim  <- run_batch(dgp_gen(FALSE), grid,
                      num_sim   = 1000,
                      base_seed = 204,
                      out_file  = "results/dimensionDGP2.csv.gz",
                      meta = list(label = "DGP2dimension"))
res_dimm <- run_batch(dgp_gen(TRUE), grid,
                      num_sim   = 1000,
                      base_seed = 204,
                      out_file  = "results/dimensionmisspecDGP2.csv.gz",
                      meta = list(label = "DGP2dimensionmisspec"))