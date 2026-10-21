# p axis runs. p/n 0.05 to 4 at n = 500 and 1000. dgp2 s = 4 overlap 1 both specs seed 202
# diag files reuse the AIPW set draws
source(here::here("runs", "aipw", "setup.R"))

grid <- expand.grid(n = c(500, 1000), pn = c(0.05, 0.1, 0.5, 1, 2, 4))
grid <- data.frame(n = grid$n, p = as.integer(grid$pn * grid$n), overlap = 1)

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = cell$p, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

for (sp in c("correct", "misspec")) {
  f_main <- paste0("aipw_", sp, "_p_dgp2.csv.gz")
  f_diag <- paste0("aipw_diag_", sp, "_p_dgp2.csv.gz")
  run_arm(f_main, dgp_gen(sp == "misspec"), grid, reps(10000))
  if (run_arm(f_diag, dgp_gen(sp == "misspec"), grid, reps(10000), estimate_fine))
    check_pairs(f_diag, f_main, c("balnet0", "glmnetcv_hajek", "oracle_hajek"),
                by = c("n", "p"))
}
