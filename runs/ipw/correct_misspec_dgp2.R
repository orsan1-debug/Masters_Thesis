## CORRECT and MISSPEC SPEC DGP2 ---------------------------------------------------------
## DGP2, correct spec + MISSPEC, iid, p = 50, overlap = 1, 1000 reps.


source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000), overlap = 1)

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

res_cor <- run_batch(dgp_gen(FALSE), grid,
                     num_sim   = 1000,
                     base_seed = 202,
                     out_file  = "results/correctspecDGP2.csv.gz",
                     meta = list(label = "DGP2correctspec"))

res_mis <- run_batch(dgp_gen(TRUE), grid,
                     num_sim   = 1000,
                     base_seed = 202,
                     out_file  = "results/misspecDGP2.csv.gz",
                     meta = list(label = "DGP2misspec"))