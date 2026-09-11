## tests/smoke_simulate.R -----------------------------------------------------
## Smoke test for R/simulate.R driver mechanics. Uses stubs so it runs
## standalone (no balnet/glmnet). Run from project root: source("tests/smoke_simulate.R")

source(here::here("R", "simulate.R"))

## --- stubs mimicking the dgp_gen / estimate_all contract ---------------------
dgp_stub <- function(cell) list(x = rnorm(cell$n))

dgp_stub_flaky <- function(cell) {            # fails deterministically per stream
  x <- rnorm(cell$n)
  if (x[1] > 0.5) stop("flaky rep")
  list(x = x)
}

dgp_stub_cell2fail <- function(cell) {        # all reps of cell 2 fail
  if (cell$n == 200) stop("cell 2 down")
  dgp_stub(cell)
}

estimate_all <- function(dat, cv_curve = FALSE) {
  m <- matrix(c(mean(dat$x), median(dat$x), sd(dat$x)), ncol = 1,
              dimnames = list(c("est_mean", "est_median", "est_sd"), "tau_hat"))
  if (cv_curve) attr(m, "cv_curve") <- dat$x[1:3]
  m
}

grid     <- data.frame(n = c(100L, 200L))
out_file <- file.path(tempdir(), "smoke.csv.gz")
args     <- list(grid = grid, num_sim = 5, base_seed = 42, out_file = out_file,
                 cv_curve_reps = 2, workers = 2)

## --- Test 1: schema ----------------------------------------------------------
r1 <- do.call(simulate_grid, c(list(dgp_stub), args))
stopifnot(
  nrow(r1) == 2 * 5 * 3,
  identical(names(r1), c("n", "sim", "estimator", "tau_hat", "err")),
  all(r1$err == ""), all(is.finite(r1$tau_hat)),
  file.exists(sub("\\.csv\\.gz$", "_session.txt", out_file))
)
curves <- readRDS(sub("\\.csv\\.gz$", "_cvcurves.rds", out_file))
stopifnot(length(curves) == 2, all(lengths(curves) == 2),
          !any(vapply(unlist(curves, recursive = FALSE), is.null, logical(1))))
cat("Test 1 (schema) PASS\n")

## --- Test 2: reproducibility -------------------------------------------------
r2 <- do.call(simulate_grid, c(list(dgp_stub), args))
stopifnot(isTRUE(all.equal(as.data.frame(r1), as.data.frame(r2))))
cat("Test 2 (reproducibility) PASS\n")

## --- Test 3: resume ----------------------------------------------------------
## Cell 2 all-fail -> driver stops, cell 1 checkpoint persists; rerun resumes.
bad <- try(do.call(simulate_grid, c(list(dgp_stub_cell2fail), args)), silent = TRUE)
ck_dir <- sub("\\.csv\\.gz$", "_cells", out_file)
stopifnot(inherits(bad, "try-error"),
          file.exists(file.path(ck_dir, "cell_001.csv.gz")))
r3 <- do.call(simulate_grid, c(list(dgp_stub), args))
stopifnot(isTRUE(all.equal(as.data.frame(r1), as.data.frame(r3))),
          !dir.exists(ck_dir))
cat("Test 3 (resume, all-fail stop) PASS\n")

## --- Test 4: per-rep error path ---------------------------------------------
args4 <- modifyList(args, list(num_sim = 20,
                               out_file = file.path(tempdir(), "smoke_flaky.csv.gz")))
r4 <- do.call(simulate_grid, c(list(dgp_stub_flaky), args4))
failed <- r4$err != ""
stopifnot(any(failed), !all(failed),
          all(is.na(r4$tau_hat[failed])), all(is.finite(r4$tau_hat[!failed])),
          all(tapply(r4$err, paste(r4$n, r4$sim), function(e) length(unique(e)) == 1)))
cat("Test 4 (error path) PASS\n")

cat("All driver smoke tests PASS\n")

## --- Test 5 (real pipeline, run separately) ----------------------------------
## Swap stubs for source("R/dgp.R"); source("R/estimators.R"), use a 2-cell real
## grid with small n, and rerun Tests 1-3. Then confirm dgp_gen returns all
## outcomes from a single (X, W) draw (shared-draws check lives in dgp.R).