## runs/aipw/aipw_overlap.R ---------------------------------------------------
## Q6 (Q1, Q7, Q8): overlap axis, dgp2, both specs, n = 1000, p = 50, seed 202.
## Ascending grid: new streams for every cell, so rows do not reproduce
## the 4-level aipw_misspec_overlap run (compare within MCSE, not to 1e-10).
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))

n_sim   <- 10000   # set to 2 for a smoke run
out_dir <- here::here("results", "aipw")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grid <- data.frame(n = 1000, overlap = c(0.05, 0.15, 0.25, 0.5, 0.75, 1, 2))

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

res_cor <- simulate_grid(dgp_gen(FALSE), grid, num_sim = n_sim, base_seed = 202,
                         out_file = file.path(out_dir,
                                              "aipw_correct_overlap7_dgp2.csv.gz"))
res_mis <- simulate_grid(dgp_gen(TRUE), grid, num_sim = n_sim, base_seed = 202,
                         out_file = file.path(out_dir,
                                              "aipw_misspec_overlap7_dgp2.csv.gz"))