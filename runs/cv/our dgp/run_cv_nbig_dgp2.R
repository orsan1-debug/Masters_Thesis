# large n lasso. the floor's rmse gap to the path min grew with n in run_cv_n_dgp2.R
# does it keep growing at Study 1's n = 10 000. alpha 0.5 dropped as a rep costs about n / 1000x
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "estimators_cv_dgp2.R"))
source(here::here("R", "simulate.R"))
source(here::here("runs", "cv", "our dgp", "path_floor.R"))
source(here::here("runs", "cv", "our dgp", "cv_helpers.R"))
stopifnot(packageVersion("glmnet") >= "5.0",
          "cv.boot.balnet" %in% getNamespaceExports("balnet"))  # cv-experimental build

# mac time at N_SIM 2500 about 3 h for n 5000 and 6-7 h for n 10 000. N_SIM 1000 about 4 h
n_sim   <- as.integer(Sys.getenv("N_SIM", 2500))   # N_SIM = 5 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "cv", "our dgp",
                                            paste0("r", n_sim)))

#' the two overlap cells at one n. one file per n so the n = 5000 results land first
grid_at <- function(n)
  data.frame(n = n, overlap = c(0.25, 0.15), cell = seq_len(2))

grids <- list(n5000 = grid_at(5000), n10000 = grid_at(10000))
seeds <- c(n5000 = 310, n10000 = 311)   # own streams after the folds grid 309

dgp_gen <- function(cell)
  share_strengths(dgp2_cv(overlap = cell$overlap, n = cell$n))

for (fam in names(grids))
  run_cv_grid(dgp_gen, grids[[fam]], num_sim = n_sim, base_seed = seeds[[fam]],
              out_file = file.path(out_dir, paste0("cv_", fam, "_dgp2.csv.gz")))
