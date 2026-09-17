## Outcome noise SNR (Q9): outcome signal scale at n = 1000, p = 50, overlap 1, both specs.
## SNR = strength^2 * Var(f); at strength 1, linear = 1.75, quad1 = 5 (runs/ipw/SNR.R).
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))

n_sim   <- 10000
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))

grid <- data.frame(n = 1000, overlap = 1,
                   strength = c(0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1))

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap,
       strength = cell$strength)

simulate_grid(dgp_gen(FALSE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_strength_dgp2.csv.gz"))
simulate_grid(dgp_gen(TRUE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_misspec_strength_dgp2.csv.gz"))