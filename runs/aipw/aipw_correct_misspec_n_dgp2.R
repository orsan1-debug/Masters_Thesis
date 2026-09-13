library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))

out_dir <- here::here("results", "aipw_10k")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grid <- data.frame(n = 1000, overlap = c(0.25, 0.5, 1, 2))

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

res_cor <- simulate_grid(dgp_gen(FALSE), grid, num_sim = 10000, base_seed = 202,
                         out_file = file.path(out_dir, "aipw_correct_overlap_dgp2.csv.gz"))
res_mis <- simulate_grid(dgp_gen(TRUE), grid, num_sim = 10000, base_seed = 202,
                         out_file = file.path(out_dir, "aipw_misspec_overlap_dgp2.csv.gz"))