## Large-n axis, correct spec only: bias and RMSE to zero with n (efficiency).
## 1000 reps: at n = 800k the estimator SE is ~0.003, so MCSE(bias) ~ 1e-4.
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))

n_sim   <- 1000
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000, 100000, 200000, 400000, 800000),
                   overlap = 1)

dgp_gen <- function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = FALSE, overlap = cell$overlap)

simulate_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_n800k_dgp2.csv.gz"))
