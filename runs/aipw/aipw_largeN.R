## Large-n axis, correct spec only: bias and RMSE to zero with n (efficiency).
## 10 000 reps: at n = 400k the ATE estimators' SD is 0.004-0.007 (linear
## outcome), so MCSE(bias) ~ 4e-5 to 7e-5. n = 800k dropped: > 2 days per cell.
## Checkpoints: <out_dir>/aipw_correct_n400k_dgp2_cells/cell_NNN.csv.gz, one per
## finished cell; a rerun resumes from them and they are deleted on completion.
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))

n_sim   <- as.integer(Sys.getenv("N_SIM", "10000"))   # N_SIM=2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))
message("out_dir: ", out_dir)                          # check this is results/aipw

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000, 100000, 200000, 400000),
                   overlap = 1)

dgp_gen <- function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = FALSE, overlap = cell$overlap)

simulate_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_n400k_dgp2.csv.gz"))
o