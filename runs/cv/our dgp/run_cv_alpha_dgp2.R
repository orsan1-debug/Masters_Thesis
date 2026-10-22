# alpha x overlap x strength. 20 cells
# replaces the overlap 0.25 sweep r1000/cv_alpha_dgp2.csv.gz
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

grid <- expand.grid(alpha = c(0, 0.25, 0.5, 0.75, 1),
                    overlap = c(1, 0.5, 0.25, 0.15))
grid$min_ratio <- floor_min_ratio(grid$alpha)    # path floor for alpha < 1
grid$cell      <- seq_len(nrow(grid))

dgp_gen <- function(cell) share_strengths(dgp2_cv(overlap = cell$overlap))

run_cv_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 305,
            out_file = file.path(out_dir, "cv_alpha_x_overlap_dgp2.csv.gz"))
