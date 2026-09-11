## CORRECT N200K ---------------------------------------------------------
## DGP2, correct spec N200K, iid, p = 50, overlap = 1, 1000 reps.


source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))
n_sim   <- as.integer(Sys.getenv("N_SIM", "1000"))          # N_SIM=2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "ipw"))  # OUT_DIR=<tmp> keeps results/ untouched

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000, 100000, 200000), overlap = 1)

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

res_cor <- run_batch(dgp_gen(FALSE), grid,
                     num_sim   = n_sim,
                     base_seed = 202,
                     out_file  = file.path(out_dir, "correctspecDGP200k.csv.gz"),
                     dgp = "dgp2", script = "runs/ipw/correctspec_dgp2_n800k.R")


## CORRECT N400K ---------------------------------------------------------
grid400 <- data.frame(n = 400000, overlap = 1)
res_cor_400k <- run_batch(dgp_gen(FALSE), grid400,
                          num_sim   = n_sim,
                          base_seed = 203,
                          out_file  = file.path(out_dir, "correctspecDGP400k.csv.gz"),
                          dgp = "dgp2", script = "runs/ipw/correctspec_dgp2_n800k.R")

## CORRECT N800K ---------------------------------------------------------
grid800 <- data.frame(n = 800000, overlap = 1)
res_cor_800k <- run_batch(dgp_gen(FALSE), grid800,
                          num_sim   = n_sim,
                          base_seed = 204,
                          out_file  = file.path(out_dir, "correctspecDGP800k.csv.gz"),
                          dgp = "dgp2", script = "runs/ipw/correctspec_dgp2_n800k.R")