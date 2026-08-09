r
## OVERLAP + OVERLAP-MISSPEC DGP2 ---------------------------------------------
## DGP2, correct + misspec, iid, p = 50, overlap grid, 1000 reps.
## Rerun territory of failed dgp2mo (err=1 batch); driver now stores conditionMessage.
source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")

grid <- expand.grid(n = c(500, 1000, 5000, 10000),
                    overlap = c(0.25, 0.5, 0.75, 1, 2))

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

res_ov  <- run_batch(dgp_gen(FALSE), grid,
                     num_sim   = 1000,
                     base_seed = 203,
                     out_file  = "results/overlapDGP2.csv.gz",
                     meta = list(label = "DGP2overlap"))

res_ovm <- run_batch(dgp_gen(TRUE), grid,
                     num_sim   = 1000,
                     base_seed = 203,
                     out_file  = "results/overlapmisspecDGP2.csv.gz",
                     meta = list(label = "DGP2overlapmisspec"))