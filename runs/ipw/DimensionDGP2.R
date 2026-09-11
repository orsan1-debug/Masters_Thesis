## DIMENSIONALITY DGP2 --------------------------------------------------------
## DGP2, correct + misspec, iid, archived E5 p/n grid, overlap = 1, 1000 reps.
source(here::here("R", "packages.R"))
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))
n_sim   <- as.integer(Sys.getenv("N_SIM", "1000"))          # N_SIM=2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "ipw"))  # OUT_DIR=<tmp> keeps results/ untouched
grid   <- expand.grid(n = c(500, 1000),
                      pn_ratio = c(0.1, 0.5, 1, 2))
grid$p <- grid$n * grid$pn_ratio
dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = cell$p, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = 1)
res_dim  <- run_batch(dgp_gen(FALSE), grid,
                      num_sim   = n_sim,
                      base_seed = 204,
                      out_file  = file.path(out_dir, "dimensionDGP2.csv.gz"),
                      dgp = "dgp2", script = "runs/ipw/DimensionDGP2.R")
res_dimm <- run_batch(dgp_gen(TRUE), grid,
                      num_sim   = n_sim,
                      base_seed = 204,
                      out_file  = file.path(out_dir, "dimensionmisspecDGP2.csv.gz"),
                      dgp = "dgp2", script = "runs/ipw/DimensionDGP2.R")