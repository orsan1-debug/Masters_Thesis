## HIGH-P SPARSITY DGP2 --------------------------------------------------------
## DGP2, correct spec, iid, n = 1000, p x s grid, decay {0,1}, overlap = 1, 500 reps.
## decay 0 = flat profile (sharp lasso-omission boundary); decay 1 links to E6.
source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")
grid <- expand.grid(n = 1000,
                    p = c(50, 250, 1000, 4000),
                    s = c(4, 16, 64, 256),
                    decay = c(0, 1))
grid <- grid[grid$s <= grid$p, ]   # 26 cells
dgp_gen <- function(cell)
  dgp2(n = cell$n, p = cell$p, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1"), outcome_set = "track_s",
       decay_ps = cell$decay, decay_out = cell$decay,
       covcor = "iid", misspec = FALSE, overlap = 1)
res_sphp <- run_batch(dgp_gen, grid,
                      num_sim   = 1000,
                      base_seed = 208,
                      out_file  = "results/sparsityhighpDGP2.csv.gz",
                      meta = list(label = "DGP2sparsityhighp"))
