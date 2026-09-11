## INSTRUMENTS DGP2 -------------------------------------------------------------
## DGP2, correct spec, iid, p = 50, overlap = 1, 1000 reps.
## Confounders fixed at 4 (outcome loads on 1..4 only); instruments are PS-only
## covariates 5..s, s = 4 + s_instr. decay is the instrument-strength lever
## (0.5 flattens weights, so instruments carry a larger eta share).
## s_instr = 0 arm is the no-instrument negative control.
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "registry.R"))
grid <- expand.grid(n = c(500, 1000, 5000, 10000),
                    s_instr = c(0, 4, 16),
                    decay = c(1, 0.5))
grid$s <- 4 + grid$s_instr
dgp_gen <- function(cell)
  dgp2(n = cell$n, p = 50, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1"),
       outcome_set = "fixed4", decay_ps = cell$decay,
       covcor = "iid", misspec = FALSE, overlap = 1)
res_iv <- run_batch(dgp_gen, grid,
                    num_sim   = 1000,
                    base_seed = 207,
                    out_file  = here::here("results", "ipw", "instrumentsDGP2.csv.gz"),
                    meta = list(label = "DGP2instruments"))




