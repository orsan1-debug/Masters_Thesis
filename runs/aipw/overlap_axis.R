# overlap axis runs. dgp2 p = 50 s = 4 both specs seed 202
# sd(eta) = 1 / overlap so smaller = worse overlap
source(here::here("runs", "aipw", "setup.R"))

dgp_gen <- function(misspec) function(cell)
  dgp2(n = cell$n, p = 50, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"),
       covcor = "iid", misspec = misspec, overlap = cell$overlap)

ov7      <- c(0.05, 0.15, 0.25, 0.5, 0.75, 1, 2)
grid_1k  <- data.frame(n = 1000,  overlap = ov7)
grid_10k <- data.frame(n = 10000, overlap = ov7)

# 7 levels at n = 1000 and 10000. AIPW set
run_arm("aipw_correct_overlap7_dgp2.csv.gz", dgp_gen(FALSE), grid_1k,  reps(10000))
run_arm("aipw_misspec_overlap7_dgp2.csv.gz", dgp_gen(TRUE),  grid_1k,  reps(10000))
run_arm("aipw_correct_overlap8_dgp2.csv.gz", dgp_gen(FALSE), grid_10k, reps(10000))
run_arm("aipw_misspec_overlap8_dgp2.csv.gz", dgp_gen(TRUE),  grid_10k, reps(10000))

# diagnostics on the n = 1000 draws
if (run_arm("aipw_diag_correct_overlap7_dgp2.csv.gz", dgp_gen(FALSE), grid_1k,
            reps(10000), estimate_fine))
  check_pairs("aipw_diag_correct_overlap7_dgp2.csv.gz",
              "aipw_correct_overlap7_dgp2.csv.gz",
              c("balnet0", "glmnetcv_hajek", "oracle_hajek"),
              by = c("n", "overlap"))
if (run_arm("aipw_full_misspec_overlap7_dgp2.csv.gz", dgp_gen(TRUE), grid_1k,
            reps(10000), estimate_full))
  check_pairs("aipw_full_misspec_overlap7_dgp2.csv.gz",
              "aipw_misspec_overlap7_dgp2.csv.gz",
              c("balnet0", "glmnetcv_hajek", "oracle_hajek"),
              by = c("n", "overlap"))

# n x overlap. 12 n levels (100-5000) x 13 overlap levels (0.01-1) 1000 reps
# own streams so cells shared w/ the runs above only agree within MCSE
n_lev   <- round(signif(exp(seq(log(100), log(5000), length.out = 12)), 2))
ov_lev  <- signif(10^seq(-2, 0, by = 1 / 6), 3)
grid_nc <- expand.grid(n = n_lev, overlap = ov_lev)
run_arm("aipw_misspec_n_by_overlap_dgp2.csv.gz", dgp_gen(TRUE),  grid_nc, reps(1000))
run_arm("aipw_correct_n_by_overlap_dgp2.csv.gz", dgp_gen(FALSE), grid_nc, reps(1000))

# first overlap run. 4 levels at n = 1000 (superseded by overlap7)
grid_4 <- data.frame(n = 1000, overlap = c(0.25, 0.5, 1, 2))
run_arm("aipw_correct_overlap_dgp2.csv.gz", dgp_gen(FALSE), grid_4, reps(10000))
run_arm("aipw_misspec_overlap_dgp2.csv.gz", dgp_gen(TRUE),  grid_4, reps(10000))
