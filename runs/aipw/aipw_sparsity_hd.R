## runs/aipw/aipw_sparsity_hd.R ----------------------------------------------
## Sparsity x dimension x correlation (Q5): p/n = 0.1, 0.5, 1 at
## n = 500, 1000 (same relative design at both n); s = 4, 16, 64 active
## covariates (s <= p; Tan 2020 RCAL rate |S| log(p)/n); iid and
## ar1 (rho 0.5, Zhao 2017 s.6.2 Sigma); decay 1, 0.5; correct spec only
## (misspec requires s = 4); linear + quad1, track_s so every active PS
## covariate is a confounder (fixed4 at s > 4 would add instruments), overlap 1.
## Mirrors runs/ipw/Sparsity.R and sparsity_2.R with the AIPW set, glmnet 5.
## p/n = 1 is the Tan 2020 AoS design point (p = 1000, n = 800).
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
stopifnot(packageVersion("glmnet") >= "5.0")

n_sim   <- as.integer(Sys.getenv("N_SIM", 1000))   # N_SIM = 2 for a smoke run
out_dir <- Sys.getenv("OUT_DIR", here::here("results", "aipw"))

grid <- expand.grid(covcor = c("iid", "ar1"), decay = c(1, 0.5),
                    s = c(4, 16, 64), pn = c(0.1, 0.5, 1),
                    n = c(500, 1000), stringsAsFactors = FALSE)
grid$p <- as.integer(grid$pn * grid$n)
grid <- grid[grid$s <= grid$p, c("n", "pn", "p", "s", "decay", "covcor")]
rownames(grid) <- NULL
stopifnot(nrow(grid) == 68)                        # 17 (p, s) pairs x 4

dgp_gen <- function(cell)
  dgp2(n = cell$n, p = cell$p, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1"), outcome_set = "track_s",
       decay_ps = cell$decay, decay_out = cell$decay,
       covcor = cell$covcor, misspec = FALSE, overlap = 1)

simulate_grid(dgp_gen, grid, num_sim = n_sim, base_seed = 202,
              out_file = file.path(out_dir, "aipw_correct_sparsity_hd_dgp2.csv.gz"))