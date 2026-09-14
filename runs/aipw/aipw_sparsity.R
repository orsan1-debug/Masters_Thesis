## Sparsity axis (Q5): s active covariates at n = 1000, p = 50, overlap 1;
## correct spec only (misspec = TRUE requires s == 4); outcome tracks s.
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))

n_sim   <- 10000
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))

grid <- data.frame(n = 1000, s = c(4, 8, 16, 32, 50))

dgp_gen <- function(cell)
  dgp2(n = cell$n, p = 50, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1", "exp"), outcome_set = "track_s",
       covcor = "iid", misspec = FALSE, overlap = 1)

simulate_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_sparsity_dgp2.csv.gz"))
