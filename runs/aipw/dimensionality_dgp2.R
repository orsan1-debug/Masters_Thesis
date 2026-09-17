## p axis (Q5): p/n at n = 500 and 1000, overlap 1, both specs. Old IPW
## grid (runs/ipw/DimensionDGP2.R) cut at p/n = 4, AIPW set, glmnet 5, 10k reps.
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
stopifnot(packageVersion("glmnet") >= "5.0")

n_sim   <- 10000
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))

grid <- expand.grid(n = c(500, 1000), pn = c(0.05, 0.1, 0.5, 1, 2, 3))
grid <- data.frame(n = grid$n, p = as.integer(grid$pn * grid$n), overlap = 1)

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = cell$p, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

simulate_grid(dgp_gen(FALSE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_p_dgp2.csv.gz"))
simulate_grid(dgp_gen(TRUE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_misspec_p_dgp2.csv.gz"))