# misspec twin of the overlap and strength grids. dgp2 misspec = TRUE (Kang-Schafer transforms)
# where is the floor not rmse optimal and which selector wins. out to overlap 0.05 strength 0
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

grid <- data.frame(overlap = c(1, 0.75, 0.5, 0.25, 0.15, 0.05), cell = 1:6)

dgp_gen <- function(cell)
  share_strengths(dgp2_cv(overlap = cell$overlap, misspec = TRUE))

run_cv_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 308,
            out_file = file.path(out_dir, "cv_misspec_dgp2.csv.gz"))
