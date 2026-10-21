# treated share runs. P(W = 1) 0.05 to 0.5 at n = 1000 and 10000. HT vs Hajek w/ rare treated arm
# dgp2 p = 50 s = 4 overlap 1 both specs seed 202
source(here::here("runs", "aipw", "setup.R"))

grid <- expand.grid(treat_prop = c(0.05, 0.1, 0.2, 0.5), n = c(1000, 10000))
grid$overlap <- 1

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap,
       treat_prop = cell$treat_prop)

run_arm("aipw_correct_treat_prop_dgp2.csv.gz", dgp_gen(FALSE), grid, reps(10000))
run_arm("aipw_misspec_treat_prop_dgp2.csv.gz", dgp_gen(TRUE),  grid, reps(10000))
