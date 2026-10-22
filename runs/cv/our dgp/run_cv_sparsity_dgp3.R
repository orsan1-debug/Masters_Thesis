# dense confounding x correlation. twin of the Study 1 sparsity grid runs/aipw/sparsity_axis.R
# other cv grids have 4 confounders in 50 iid covs. interior optimum may be a sparsity artefact
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "estimators_cv_dgp2.R"))
source(here::here("R", "simulate.R"))
source(here::here("runs", "cv", "our dgp", "path_floor.R"))
source(here::here("runs", "cv", "our dgp", "cv_helpers.R"))
stopifnot(packageVersion("glmnet") >= "5.0",
          "cv.boot.balnet" %in% getNamespaceExports("balnet"),  # in balnet 0.0.4 from pkg/
          exists("dgp3"))

# mac time at N_SIM 2500 about 4-6 h
n_sim   <- as.integer(Sys.getenv("N_SIM", 5))   # N_SIM = 5 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "cv", "our dgp",
                                            paste0("r", n_sim)))

grid <- expand.grid(s = c(4, 16, 64), covcor = c("iid", "ar1", "equi"),
                    stringsAsFactors = FALSE)
grid$n       <- 1000
grid$p       <- 100
grid$decay   <- 0.5
grid$overlap <- 0.25   # where the n grid's interior optima are clear
grid$cell    <- seq_len(nrow(grid))   # adding cells later would move these streams

#' dgp_correct() of runs/aipw/sparsity_axis.R via dgp3() for covcor equi. plus the exp outcome
#' iid or ar1 with misspec FALSE gives the same draws as dgp2()
dgp3_cv <- function(cell)
  dgp3(n = cell$n, p = cell$p, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1", "exp"), outcome_set = "track_s",
       decay_ps = cell$decay, decay_out = cell$decay,
       covcor = cell$covcor, misspec = FALSE, overlap = cell$overlap)

dgp_gen <- function(cell) share_strengths(dgp3_cv(cell))

# cv_sparsity_misspec_dgp3 was this run w/ misspec = TRUE and misspec_form tan_obs in dgp3_cv()
run_cv_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 312,
            out_file = file.path(out_dir, "cv_sparsity_dgp3.csv.gz"))
