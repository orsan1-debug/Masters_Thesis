# outcome signal strength runs. dgp2 n = 1000 p = 50 s = 4 overlap 1 both specs seed 202
# SNR = strength^2 * Var(f). at strength 1 linear = 1.75 and quad1 = 5
source(here::here("runs", "aipw", "setup.R"))

grid <- data.frame(n = 1000, overlap = 1,
                   strength = c(0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1))

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap,
       strength = cell$strength)

run_arm("aipw_correct_strength_dgp2.csv.gz", dgp_gen(FALSE), grid, reps(10000))
run_arm("aipw_misspec_strength_dgp2.csv.gz", dgp_gen(TRUE),  grid, reps(10000))
