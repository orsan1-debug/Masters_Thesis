# Wong and Chan check (App H). n = 5000. noise sd 1 10 30 100 x logit multiplier 1 to 4. 1000 reps
#   Rscript replication/runs/wong_chan.R   16 cells to results/cv/wz_replication. ~6 s a rep a core

source(here::here("replication", "R", "wong_chan.R"))
source(here::here("replication", "R", "simulate.R"))

# thesis matches the multiplier 1 2 and 4 cells to its misspec grid
run_wong_chan_cells(
  cells        = local({ g <- expand.grid(sigma = c(1, 10, 30, 100), c = c(1, 2, 3, 4)); g[cells_to_run(g), ] }),
  draw         = function(cell) generate_wong_chan(5000, overlap = cell$c, sigma = cell$sigma),
  tolerances   = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1),
  replications = min(as.integer(Sys.getenv("N_SIM", 1000)), 1000L),   # no more than the stored 1000
  out_dir      = Sys.getenv("OUT_DIR", here::here("results", "cv", "wz_replication")),
  file_of      = function(cell) sprintf("wc_noise_s%g_c%g.rds", cell$sigma, cell$c),
  workers      = as.integer(Sys.getenv("WORKERS", max(1L, parallel::detectCores() - 1L))))
