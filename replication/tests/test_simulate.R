# seed streams and run_grid() vs the record (R/simulate.R) and the stored misspec sample size draws
# also block columns and specs sharing a block. w/o results/ stored draws come from tests/fixtures/

library(testthat)
source(here::here("replication", "R", "dgp.R"))
source(here::here("replication", "R", "estimators.R"))
source(here::here("replication", "R", "simulate.R"))

record <- new.env()
sys.source(here::here("R", "dgp.R"), envir = record)
sys.source(here::here("R", "estimators_ipw.R"), envir = record)
sys.source(here::here("R", "simulate.R"), envir = record)

#' a stored file or if results/ is absent its rows in tests/fixtures/ (from make_fixtures.R)
stored_path <- function(...) {
  path <- here::here("results", ...)
  ifelse(file.exists(path), path, here::here("replication", "tests", "fixtures", ...))
}

test_that("the seed streams are the record's", {
  expect_identical(seed_streams(202, 25), record$make_seeds(202, 25))
  expect_identical(seed_streams(301, 7), record$make_seeds(301, 7))
})

test_that("stored draws are regenerated: misspecified sample-size run, n = 500, first three replications", {
  f <- stored_path("aipw", "aipw_misspec_n_dgp2.csv.gz")
  if (!file.exists(f)) skip("results file not present")
  stored <- data.table::fread(f)[n == 500 & sim <= 3]
  streams <- seed_streams(202, 3)                   # cell 1 uses streams 1 to 3
  for (r in 1:3) {
    assign(".Random.seed", streams[[r]], envir = globalenv())
    d <- generate_sample(n = 500, p = 50, misspecification = "kang_schafer")
    ref <- estimate_reference(d)
    bal <- fit_balancing(d$X, d$W)
    balnet0 <- colMeans(d$Y * drop(bal$weights$treated - bal$weights$control))
    for (est in c("oracle_hajek", "balnet0")) {
      got <- if (est == "oracle_hajek") ref["oracle_hajek", ] else balnet0
      want <- unlist(stored[sim == r & estimator == est, .(linear, quad1, exp)])
      expect_equal(unname(got), unname(want), tolerance = 1e-10, label = paste(est, "rep", r))
    }
    # MLE rows can depend on the glmnet version so reported not asserted
    mle <- fit_mle(d$X, d$W)
    got <- ate_hajek(d$Y, d$W, mle$e_hat, 1 - mle$e_hat)
    want <- unlist(stored[sim == r & estimator == "glmnetcv_hajek", .(linear, quad1, exp)])
    cat(sprintf("  rep %d glmnetcv_hajek, largest difference to stored: %.2e (glmnet %s here)\n",
                r, max(abs(got - want)), packageVersion("glmnet")))
  }
})

test_that("run_grid() reproduces the record's driver on a small grid", {
  grid <- data.frame(n = c(200, 300), overlap = c(1, 0.5))
  gen_new <- function(cell) generate_sample(n = cell$n, p = 10, c = cell$overlap)
  gen_old <- function(cell) record$dgp2(n = cell$n, p = 10, s = 4, signs = "pos",
                                        outcome = c("linear", "quad1", "exp"),
                                        covcor = "iid", overlap = cell$overlap)
  tmp <- tempfile(); dir.create(tmp)
  new <- run_grid(gen_new, grid, replications = 3, base_seed = 7,
                  out_file = file.path(tmp, "new.csv.gz"), estimate = estimate_reference)
  old_opts <- options(sim.plan = future::sequential)
  on.exit(options(old_opts), add = TRUE)
  old <- suppressWarnings(
    record$simulate_grid(gen_old, grid, num_sim = 3, base_seed = 7,
                         out_file = file.path(tmp, "old.csv.gz"),
                         estimator = record$estimate_add, workers = 1))
  data.table::setnames(old, c("quad1", "exp"), c("quadratic", "exponential"))
  expect_equal(as.data.frame(new), as.data.frame(old), tolerance = 1e-12)
  expect_true(file.exists(file.path(tmp, "new_session.txt")))
})

test_that("a block column gives a cell the streams of the grid it was first run in", {
  draw <- function(cell) generate_sample(n = cell$n, p = 10)
  tmp <- tempfile(); dir.create(tmp)
  run <- function(grid, file) as.data.frame(run_grid(draw, grid, replications = 2, base_seed = 7,
    out_file = file.path(tmp, file), estimate = estimate_reference, stream_block = 5))
  first  <- run(data.frame(n = c(200, 300)), "first.csv.gz")
  second <- run(data.frame(n = 400), "second.csv.gz")
  merged <- run(data.frame(n = c(200, 300, 400), block = c(1, 2, 1)), "merged.csv.gz")
  expect_equal(merged, rbind(first, second), tolerance = 0)
})

test_that("two specifications in one block give the numbers of two separate runs", {
  draw <- function(cell)
    generate_sample(n = 300, p = 10, misspecification = if (cell$spec == "misspecified") "kang_schafer" else "none")
  tmp <- tempfile(); dir.create(tmp)
  run <- function(grid, file) {
    out <- as.data.frame(run_grid(draw, grid, replications = 2, base_seed = 202, out_file = file.path(tmp, file),
                                  estimate = estimate_study1_with_reference))
    out <- out[!grepl("^time_", out$estimator), ]
    rownames(out) <- NULL
    out
  }
  both  <- run(data.frame(spec = c("correct", "misspecified"), block = 1), "both.csv.gz")
  alone <- rbind(run(data.frame(spec = "correct"), "correct.csv.gz"),
                 run(data.frame(spec = "misspecified"), "misspecified.csv.gz"))
  expect_equal(both, alone, tolerance = 0)
})
