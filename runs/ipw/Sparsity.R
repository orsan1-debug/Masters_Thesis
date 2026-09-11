## SPARSITY DGP2 ---------------------------------------------------------------
## DGP2, correct spec, iid, p = 50, s x decay grid, overlap = 1, 1000 reps.
## s = 50 = p is the exactly dense case; decay 0.5 flattens the coefficient
## profile; track_s, outcome shares weights. Supersedes failed dgp2spars batch.
source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")
grid <- expand.grid(n = c(500, 1000, 5000, 10000),
                    s = c(4, 8, 16, 32, 50),
                    decay = c(1, 0.5),
                    covcor = c("iid", "ar1"),
                    stringsAsFactors = FALSE)
dgp_gen <- function(cell)
  dgp2(n = cell$n, p = 50, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1"),
       outcome_set = "track_s", decay_ps = cell$decay, decay_out = cell$decay,
       covcor = cell$covcor, misspec = FALSE, overlap = 1)
res_sp <- run_batch(dgp_gen, grid,
                    num_sim   = 1000,
                    base_seed = 206,
                    out_file  = "results/sparsityDGP2.csv.gz",
                    meta = list(label = "DGP2sparsity"))




