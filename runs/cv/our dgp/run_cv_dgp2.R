# cv twin of the basic_dgp grids. lasso correct spec. overlap axis strength axis and both crossed
# each strength is its own draw here (the other cv runs share one draw)
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

grids <- list(
  overlap         = data.frame(overlap = c(1, 0.75, 0.5, 0.25, 0.15, 0.05),
                               strength = 1),
  noise           = data.frame(overlap = 1, strength = cv_strengths),
  noise_x_overlap = expand.grid(overlap = c(1, 0.5, 0.25, 0.15),
                                strength = cv_strengths))
seeds <- c(overlap = 301, noise = 302, noise_x_overlap = 303)

dgp_gen <- function(cell)
  dgp2_cv(overlap = cell$overlap, strength = cell$strength)

for (fam in names(grids))
  run_cv_grid(dgp_gen, cbind(grids[[fam]], cell = seq_len(nrow(grids[[fam]]))),
              num_sim = n_sim, base_seed = seeds[[fam]],
              out_file = file.path(out_dir, paste0("cv_", fam, "_dgp2.csv.gz")))
