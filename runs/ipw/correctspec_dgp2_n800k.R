## CORRECT N200K ---------------------------------------------------------
## DGP2, correct spec N200K, iid, p = 50, overlap = 1, 1000 reps.


source("R/dgp.R")
source("R/estimators.R")
source("R/simulate.R")
source("R/registry.R")

grid <- data.frame(n = c(500, 1000, 5000, 10000, 50000, 100000, 200000), overlap = 1)

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

res_cor <- run_batch(dgp_gen(FALSE), grid,
                     num_sim   = 1000,
                     base_seed = 202,
                     out_file  = "results/correctspecDGP200k.csv.gz",
                     meta = list(label = "DGP200kcorrectspec"))


## CORRECT N400K ---------------------------------------------------------
grid400 <- data.frame(n = 400000, overlap = 1)
res_cor_400k <- run_batch(dgp_gen(FALSE), grid400,
                          num_sim   = 1000,
                          base_seed = 203,
                          out_file  = "results/correctspecDGP400k.csv.gz",
                          meta = list(label = "DGP400kcorrectspec"))

## CORRECT N800K ---------------------------------------------------------
grid800 <- data.frame(n = 800000, overlap = 1)
res_cor_800k <- run_batch(dgp_gen(FALSE), grid800,
                          num_sim   = 1000,
                          base_seed = 204,
                          out_file  = "results/correctspecDGP800k.csv.gz",
                          meta = list(label = "DGP800kcorrectspec"))