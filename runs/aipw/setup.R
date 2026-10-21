# shared setup for runs/aipw. one script per design axis and one results file per run_arm()
# all runs seed 202. stream = (cell - 1) * num_sim + rep so same grid order and num_sim pair by sim
library(balnet)
library(glmnet)
source(here::here("R", "dgp.R"))
source(here::here("R", "estimators_ipw.R"))
source(here::here("R", "simulate.R"))
stopifnot(packageVersion("glmnet") >= "5.0")             # same as the stored runs
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")                  # 1 BLAS thread per worker on mac

# smoke run = N_SIM=2 OUT_DIR=<temp>. only cell 1 then pairs with the stored files
res_dir <- here::here("results", "aipw")
out_dir <- Sys.getenv("OUT_DIR", res_dir)
smoke   <- out_dir != res_dir
n_workers <- as.integer(Sys.getenv(
  "WORKERS", max(1L, parallelly::availableCores(logical = FALSE) - 1L)))
message("out_dir: ", out_dir)

#' reps per cell or N_SIM on a smoke run
reps <- function(n) as.integer(Sys.getenv("N_SIM", n))

#' run one arm to one results file. skips if the file exists or the name doesnt match ARMS
#' grid row order sets the seed streams so never reorder or extend a stored grid
run_arm <- function(file, dgp_gen, grid, num_sim, estimator = estimate_all,
                    workers = n_workers) {
  out <- file.path(out_dir, file)
  if (file.exists(out) || !grepl(Sys.getenv("ARMS", "."), file)) {
    message("skip ", file)
    return(invisible(FALSE))
  }
  simulate_grid(dgp_gen, grid, num_sim = num_sim, base_seed = 202,
                out_file = out, estimator = estimator, workers = workers)
  invisible(TRUE)
}

#' check two results files hold the same draws. stops on a full run if max diff > tol
#' draw only rows (oracle_hajek unadjusted oracle_aipw) agree to 1e-12 and refits to ~1e-6
check_pairs <- function(a, b, arms, by = "n", tol = 1e-6) {
  read <- function(f) {
    path <- file.path(out_dir, f)
    if (!file.exists(path)) path <- file.path(res_dir, f)
    data.table::fread(path, colClasses = list(character = "err"))[
      estimator %in% arms]
  }
  x <- read(a)
  y <- read(b)
  outcomes <- intersect(c("linear", "quad1", "exp"), names(x))
  m <- merge(x, y, by = c(by, "sim", "estimator"))
  m$diff <- do.call(pmax, lapply(outcomes, function(o)
    abs(m[[paste0(o, ".x")]] - m[[paste0(o, ".y")]])))
  chk <- m[, .(reps = .N, max_diff = max(diff, na.rm = TRUE)),
           by = c(by, "estimator")]
  cat(a, "vs", b, "\n")
  print(chk)
  if (!smoke) stopifnot(chk$max_diff < tol)
  invisible(chk)
}

#' add-only rows then diag ps fits. estimate_add() has no RNG so shared rows pair w/ AIPW files
#' ols is NA when an arm has fewer than p + 1 units (lm.fit rank deficient)
estimate_fine <- function(data) {
  add <- estimate_add(data)[c("ols", "unadjusted", "oracle_aipw"), ]
  if (ncol(data$X) + 1 > min(table(data$W))) add["ols", ] <- NA
  rbind(add, estimate_diag(data))
}

#' add-only rows + AIPW set + diagnostics in one pass
#' diag restarts from estimate_all()'s RNG state so same cv folds and shared rows must match
estimate_full <- function(data) {
  force(data)                                      # draw runs here
  seed  <- get(".Random.seed", envir = globalenv())
  a_all <- estimate_all(data)
  assign(".Random.seed", seed, envir = globalenv())
  a_diag <- estimate_diag(data)
  shared <- setdiff(intersect(rownames(a_all), rownames(a_diag)),
                    c("time_bal", "time_ps"))
  stopifnot(identical(a_all[shared, ], a_diag[shared, ]))
  rbind(estimate_add(data)[c("ols", "unadjusted", "oracle_aipw"), ], a_all,
        a_diag[setdiff(rownames(a_diag), rownames(a_all)), ])
}
