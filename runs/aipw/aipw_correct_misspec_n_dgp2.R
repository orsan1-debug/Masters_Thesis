## n axis (Q2): n at p = 50, overlap 1, both specs. Old IPW grid
## (runs/ipw/correct_misspec_dgp2.R), AIPW set, glmnet 5, 10k reps.
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
stopifnot(packageVersion("glmnet") >= "5.0")

n_sim   <- 10000
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000), overlap = 1)

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

simulate_grid(dgp_gen(FALSE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_n_dgp2.csv.gz"))
simulate_grid(dgp_gen(TRUE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_misspec_n_dgp2.csv.gz"))