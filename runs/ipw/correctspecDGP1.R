## CORRECT SPEC DGP1 ---------------------------------------------------------
## DGP1, correct spec, iid, p = 50, overlap = 1, 1000 reps.
## Supersedes legacy E1 + REDUX. Floor 1e-4 via estimate_all defaults.
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000), overlap = 1)

dgp_gen_e1 <- function(cell)
  dgp1(n = cell$n, p = 50, outcome = c("linear", "quad1", "exp"),
       overlap = cell$overlap)

res <- run_batch(dgp_gen_e1, grid,
                 num_sim   = 1000,
                 base_seed = 101,
                 out_file  = here::here("results", "ipw", "correctspecDGP1.csv.gz"),
                 meta = list(label = "DGP1correctspec"))
