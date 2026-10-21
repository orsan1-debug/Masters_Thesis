# p x overlap x penalty (twin of basic_dgp dima/dimhi). interior optimum as p nears arm size?
# in basic_dgp selectors won by 7-23% at p = 500 (R^2 0.02) and the gap closed at p = 1000
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

n_sim   <- as.integer(Sys.getenv("N_SIM", 2500))   # N_SIM = 5 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "cv", "our dgp",
                                            paste0("r", n_sim)))

grid <- expand.grid(p = c(50, 250, 500, 1000), overlap = c(0.5, 0.25),
                    alpha = c(0.5, 1))
grid$min_ratio <- floor_min_ratio(grid$alpha)    # path floor for alpha < 1
grid$cell      <- seq_len(nrow(grid))

dgp_gen <- function(cell)
  share_strengths(dgp2_cv(overlap = cell$overlap, p = cell$p))

run_cv_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 307,
            out_file = file.path(out_dir, "cv_p_dgp2.csv.gz"))
