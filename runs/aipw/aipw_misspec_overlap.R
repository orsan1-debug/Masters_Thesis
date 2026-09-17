library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))

n_sim   <- as.integer(Sys.getenv("N_SIM", 10000))   # N_SIM = 2 for a smoke
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))

grid <- data.frame(n = 1000, overlap = c(0.05, 0.15, 0.25, 0.5, 0.75, 1, 2))

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

simulate_grid(dgp_gen(FALSE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_overlap_dgp2.csv.gz"))
simulate_grid(dgp_gen(TRUE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_misspec_overlap_dgp2.csv.gz"))
