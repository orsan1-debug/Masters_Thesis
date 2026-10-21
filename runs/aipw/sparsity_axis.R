# sparsity x dim x correlation runs. s = 4 16 64 (s <= p) p/n 0.1 0.5 1 at n = 500 and 1000
# decay 1 and 0.5. overlap 1 1000 reps seed 202. p/n = 1 is the Tan 2020 AoS design point
source(here::here("runs", "aipw", "setup.R"))
stopifnot(exists("dgp3"))

#' sparsity grid for the given covcor. iid or ar1 (rho 0.5) or equi (Sigma_jk = 0.5)
grid_hd <- function(covcor) {
  g <- expand.grid(covcor = covcor, decay = c(1, 0.5),
                   s = c(4, 16, 64), pn = c(0.1, 0.5, 1),
                   n = c(500, 1000), stringsAsFactors = FALSE)
  g$p <- as.integer(g$pn * g$n)
  g <- g[g$s <= g$p, c("n", "pn", "p", "s", "decay", "covcor")]
  rownames(g) <- NULL
  g
}

# linear and quad1. track_s = every active ps covariate is a confounder
dgp_correct <- function(cell)
  dgp2(n = cell$n, p = cell$p, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1"), outcome_set = "track_s",
       decay_ps = cell$decay, decay_out = cell$decay,
       covcor = cell$covcor, misspec = FALSE, overlap = 1)

# misspec is dgp3 tan_obs. x_dag = X + {(X + 1)_+}^2 on all s active cols (Tan 2020 AoS s4.2)
# goes past the s = 4 limit of Kang-Schafer. truth W Y are dgp2's so pairs rep by rep w/ correct
dgp_tan <- function(misspec) function(cell)
  dgp3(n = cell$n, p = cell$p, s = cell$s, signs = "pos",
       outcome = c("linear", "quad1"), outcome_set = "track_s",
       decay_ps = cell$decay, decay_out = cell$decay,
       covcor = cell$covcor, misspec = misspec, misspec_form = "tan_obs",
       overlap = 1)

# iid and ar1
grid <- grid_hd(c("iid", "ar1"))
stopifnot(nrow(grid) == 68)
run_arm("aipw_correct_sparsity_hd_dgp2.csv.gz", dgp_correct, grid, reps(1000))
run_arm("aipw_misspec_tanobs_sparsity_hd_dgp3.csv.gz", dgp_tan(TRUE), grid,
        reps(1000))

# equicorrelated on its own grid. adding equi above would reorder its cells and seed streams
grid_equi <- grid_hd("equi")
stopifnot(nrow(grid_equi) == 34)
run_arm("aipw_correct_sparsity_hd_equi_dgp3.csv.gz", dgp_tan(FALSE), grid_equi,
        reps(1000))
run_arm("aipw_misspec_tanobs_sparsity_hd_equi_dgp3.csv.gz", dgp_tan(TRUE),
        grid_equi, reps(1000))
