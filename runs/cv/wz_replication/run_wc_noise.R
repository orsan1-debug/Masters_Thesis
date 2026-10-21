# noise sd x overlap on the Wong & Chan design. ATT n = 5000. sigma 1 c 1 is the paper's design
# interior best lambda on model A as noise grows? Ben-Michael 2021 s7.3 Bruns-Smith 2025 s7.1.1
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_cv.R"))
source(here::here("R", "registry.R"))
source(here::here("runs", "cv", "wz_replication", "path_floor.R"))
source(here::here("runs", "cv", "wz_replication", "wc_helpers.R"))

n_rep   <- as.integer(Sys.getenv("N_SIM", 1000))   # N_SIM = 2 for a smoke run
res_dir <- Sys.getenv("OUT_DIR", here::here("results", "cv", "wz_replication"))

cells <- expand.grid(sigma = c(1, 10, 30, 100), c = c(1, 2, 3, 4))
grid  <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)  # Algorithm 1

draw <- function(cell) dgp_wc_overlap(5000, overlap = cell$c, sigma = cell$sigma)

run_wc_cells(cells, draw, grid = grid, n_rep = n_rep, res_dir = res_dir,
             file_of = function(cell)
               paste0("wc_noise_s", cell$sigma, "_c", cell$c, ".rds"),
             glob = "wc_noise_s*_c*.rds",
             script = "runs/cv/wz_replication/run_wc_noise.R")
