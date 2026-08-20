## TREAT PROP GRID ---------------------------------------------------------------
source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")
grid <- expand.grid(n = c(1000, 10000), overlap = 1, strength = 1,
                    treat_prop = c(0.05, 0.1, 0.2, 0.5))
dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap,
       strength = cell$strength, treat_prop = cell$treat_prop)
res_prev <- run_batch(dgp_gen(FALSE), grid,
                      num_sim = 1000, base_seed = 204,
                      out_file = "results/prevDGP2correct.csv.gz",
                      meta = list(label = "DGP2prevcorrectspec"))