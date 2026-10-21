# n axis runs. dgp2 p = 50 s = 4 overlap 1 seed 202
# grid N = first 5 cells of grid A so same draws. fine grid has its own streams
source(here::here("runs", "aipw", "setup.R"))

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

grid_n <- data.frame(n = c(500, 1000, 5000, 10000, 50000), overlap = 1)
grid_a <- data.frame(n = c(500, 1000, 5000, 10000, 50000, 100000, 200000, 400000),
                     overlap = 1)
grid_fine <- data.frame(n = c(250, 500, 1000, seq(2000, 20000, by = 1000),
                              30000, 40000, 50000), overlap = 1)    # 25 levels

# grid N both specs
run_arm("aipw_correct_n_dgp2.csv.gz", dgp_gen(FALSE), grid_n, reps(10000))
run_arm("aipw_misspec_n_dgp2.csv.gz", dgp_gen(TRUE),  grid_n, reps(10000))

# ols unadjusted oracle_aipw. oracle_hajek kept to pair with the main file
if (run_arm("aipw_add_misspec_n_dgp2.csv.gz", dgp_gen(TRUE), grid_n,
            reps(10000), estimate_add))
  check_pairs("aipw_add_misspec_n_dgp2.csv.gz", "aipw_misspec_n_dgp2.csv.gz",
              "oracle_hajek", tol = 1e-12)

# path endpoint grid. max.imbalance 1e-2 to 1e-7 by half decades
# 1e-4 not here as its the main file
for (sp in c("correct", "misspec")) {
  f <- paste0("aipw_tol_", sp, "_n_dgp2.csv.gz")
  if (run_arm(f, dgp_gen(sp == "misspec"), grid_n, reps(10000), estimate_tol))
    check_pairs(f, paste0("aipw_", sp, "_n_dgp2.csv.gz"), "oracle_hajek",
                tol = 1e-12)
}

# grid A large n correct only. 800k dropped took over 2 days a cell
run_arm("aipw_correct_n400k_dgp2.csv.gz", dgp_gen(FALSE), grid_a, reps(10000))

# max 6 workers or the 24 GB mac swaps at n >= 100k
# stored file only has oracle_aipw to 50k (added later) a rerun has all n
if (run_arm("aipw_add_correct_n400k_dgp2.csv.gz", dgp_gen(FALSE), grid_a,
            reps(10000), estimate_add, workers = min(6L, n_workers))) {
  check_pairs("aipw_add_correct_n400k_dgp2.csv.gz",
              "aipw_correct_n400k_dgp2.csv.gz", "oracle_hajek", tol = 1e-12)
  check_pairs("aipw_add_correct_n400k_dgp2.csv.gz",
              "aipw_correct_n_dgp2.csv.gz", "oracle_hajek", tol = 1e-12)
}

# n = 1 million misspec. 10 batches of 1000 reps on own streams
run_arm("aipw_misspec_n_dgp2n1mill.csv.gz", dgp_gen(TRUE),
        data.frame(n = 1e6, overlap = 1, batch = 1:10), reps(1000))

# fine grid diagnostics only (no outcome models) both specs
# specs share the draw so the oracle and unadjusted rows match
ran <- c(run_arm("aipw_diag_correct_nfine_dgp2.csv.gz", dgp_gen(FALSE),
                 grid_fine, reps(10000), estimate_fine),
         run_arm("aipw_diag_misspec_nfine_dgp2.csv.gz", dgp_gen(TRUE),
                 grid_fine, reps(10000), estimate_fine))
if (any(ran))
  check_pairs("aipw_diag_correct_nfine_dgp2.csv.gz",
              "aipw_diag_misspec_nfine_dgp2.csv.gz",
              c("oracle_hajek", "unadjusted", "oracle_aipw"), tol = 1e-12)

# misspec to 200k all estimators in one go. cells 1-25 are grid_fine
run_arm("aipw_full_misspec_nfine_dgp2.csv.gz", dgp_gen(TRUE),
        data.frame(n = c(grid_fine$n, 75000, 1e5, 1.5e5, 2e5), overlap = 1),
        reps(10000), estimate_full)

# n = 250 both specs. reuses cell 1 of the fine grid
for (sp in c("correct", "misspec")) {
  f <- paste0("aipw_", sp, "_n250_dgp2.csv.gz")
  if (run_arm(f, dgp_gen(sp == "misspec"), data.frame(n = 250, overlap = 1),
              reps(10000)))
    check_pairs(f, paste0("aipw_diag_", sp, "_nfine_dgp2.csv.gz"),
                c("balnet0", "glmnetcv_hajek", "oracle_hajek"))
}
