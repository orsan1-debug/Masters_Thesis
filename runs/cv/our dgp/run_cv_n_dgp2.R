# n x overlap x penalty (enet 0.5). 16 cells
# does selector gain over the floor grow with n at poor overlap
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

grid <- expand.grid(n = c(250, 500, 1000, 2000), overlap = c(0.25, 0.15),
                    alpha = c(0.5, 1))
grid$min_ratio <- floor_min_ratio(grid$alpha)    # path floor for alpha < 1
grid$cell      <- seq_len(nrow(grid))

dgp_gen <- function(cell)
  share_strengths(dgp2_cv(overlap = cell$overlap, n = cell$n))

run_cv_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 306,
            out_file = file.path(out_dir, "cv_n_dgp2.csv.gz"))
