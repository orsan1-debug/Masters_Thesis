# cv folds K = 3 5 10. 3 is the cv.balnet min. 10 as in the other cv runs
# same seed for every K. fold shuffle takes the same rng draws so only the fold split changes
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "estimators_cv_dgp2.R"))
source(here::here("R", "simulate.R"))
source(here::here("runs", "cv", "our dgp", "path_floor.R"))
source(here::here("runs", "cv", "our dgp", "cv_helpers.R"))
stopifnot(packageVersion("glmnet") >= "5.0",
          "cv.boot.balnet" %in% getNamespaceExports("balnet"))  # in balnet 0.0.4 from pkg/

n_sim   <- as.integer(Sys.getenv("N_SIM", 2500))   # N_SIM = 5 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "cv", "our dgp",
                                            paste0("r", n_sim)))

grid <- data.frame(overlap = c(1, 0.5, 0.25, 0.15), cell = 1:4)

dgp_gen <- function(cell) share_strengths(dgp2_cv(overlap = cell$overlap))

for (k in c(10, 5, 3))
  run_cv_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 309, nfolds = k,
              out_file = file.path(out_dir, sprintf("cv_folds_k%d_dgp2.csv.gz", k)))
