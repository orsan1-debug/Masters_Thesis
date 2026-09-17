## runs/aipw/aipw_treat_prop.R -----------------------------------------------
## Treated share axis: P(W = 1) = 0.05 ... 0.5 at n = 1000, 10000, p = 50,
## s = 4, overlap 1, both specs, AIPW set, glmnet 5. Old IPW grid
## (runs/ipw/treat_prop.R, correct only, seed 204); exp outcome added to match
## the other p = 50 AIPW batches (HT vs Hajek under a rare treated arm).
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
stopifnot(packageVersion("glmnet") >= "5.0")

n_sim   <- as.integer(Sys.getenv("N_SIM", 10000))  # N_SIM = 2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))

grid <- expand.grid(treat_prop = c(0.05, 0.1, 0.2, 0.5), n = c(1000, 10000))
grid$overlap <- 1

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap,
       treat_prop = cell$treat_prop)

simulate_grid(dgp_gen(FALSE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_treat_prop_dgp2.csv.gz"))
simulate_grid(dgp_gen(TRUE), grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_misspec_treat_prop_dgp2.csv.gz"))