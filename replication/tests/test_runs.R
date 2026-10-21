# runs/study1.R and runs/study2.R vs the stored results (draws lambda grids and a few reps). ~4 min
# w/o results/ the tests read the rows they use from tests/fixtures/ (make_fixtures.R)

library(testthat)
source(here::here("replication", "R", "dgp.R"))
source(here::here("replication", "R", "estimators.R"))
source(here::here("replication", "R", "selectors.R"))
source(here::here("replication", "R", "simulate.R"))

study1 <- local({ source(here::here("replication", "runs", "study1.R"), local = TRUE); grids })
study2 <- local({ source(here::here("replication", "runs", "study2.R"), local = TRUE); grids })
key <- function(d) do.call(paste, lapply(d, function(v) if (is.numeric(v)) sprintf("%.10g", v) else v))
#' a stored file or if results/ is absent its rows in tests/fixtures/ (from make_fixtures.R)
stored_path <- function(...) {
  path <- here::here("results", ...)
  ifelse(file.exists(path), path, here::here("replication", "tests", "fixtures", ...))
}

# oracle hajek of rep 1 checks the draws w/o any fitting
test_that("every Study 1 grid draws the stored samples (first replication of every cell)", {
  for (name in names(study1)) {
    g <- study1[[name]]
    block   <- stream_blocks(g$cells)
    design  <- g$cells[setdiff(names(g$cells), "block")]
    point   <- setdiff(names(design), "spec")
    streams <- seed_streams(202, (max(block) - 1) * g$replications + 1)
    for (spec in names(g$stored)) {
      files <- stored_path("aipw", paste0(g$stored[[spec]], ".csv.gz"))
      files <- files[file.exists(files)]
      if (!length(files)) next
      rows  <- which(design$spec == spec)
      drawn <- do.call(rbind, lapply(rows, function(i) {
        assign(".Random.seed", streams[[(block[i] - 1) * g$replications + 1]], envir = globalenv())
        d <- g$draw(as.list(design[i, , drop = FALSE]))
        ate_hajek(as.matrix(d$Y), d$W, d$e, d$e0)
      }))
      rownames(drawn) <- key(design[rows, point, drop = FALSE])
      for (f in files) {
        old <- data.table::fread(f)[estimator == "oracle_hajek" & sim == 1]
        if ("batch" %in% names(old)) old <- old[batch == 1]
        cells <- key(old[, point, with = FALSE])
        expect_true(all(cells %in% rownames(drawn)), label = paste(name, basename(f), "cells in the grid"))
        expect_equal(unname(drawn[cells, , drop = FALSE]),
                     unname(as.matrix(old[, intersect(c("linear", "quad1", "exp"), names(old)), with = FALSE])),
                     tolerance = 1e-10, label = paste(name, basename(f)))
      }
    }
  }
})

test_that("study1.R pairs with the stored sample-size files, both specifications, every row", {
  stored <- c(correct      = "aipw_correct_n_dgp2",  correct      = "aipw_add_correct_n400k_dgp2",
              misspecified = "aipw_misspec_n_dgp2",  misspecified = "aipw_add_misspec_n_dgp2")
  paths <- stored_path("aipw", paste0(stored, ".csv.gz"))
  if (!all(file.exists(paths))) skip("results files not present")
  tmp <- tempfile("run"); dir.create(tmp)
  withr::local_envvar(c(N_SIM = "2", OUT_DIR = tmp, WORKERS = "1", GRIDS = "^sample_size$",
                        MAX_N = "10000", CELLS = NA))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(here::here("replication", "runs", "study1.R")),
                    stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0L)
  new <- data.table::fread(file.path(tmp, "sample_size.csv.gz"))[!grepl("^time_", estimator)]
  data.table::setnames(new, c("quadratic", "exponential"), c("quad1", "exp"))
  matched <- NULL
  for (k in seq_along(stored)) {
    old  <- data.table::fread(paths[k])[sim <= 2 & n <= 10000 & !grepl("^time_", estimator)]
    both <- merge(new[spec == names(stored)[k]], old, by = c("n", "sim", "estimator"), suffixes = c("", ".stored"))
    expect_equal(nrow(both), nrow(old), label = paste(stored[k], "rows matched"))
    expect_equal(unname(as.matrix(both[, .(linear, quad1, exp)])),
                 unname(as.matrix(both[, .(linear.stored, quad1.stored, exp.stored)])),
                 tolerance = 1e-8, label = stored[k])
    matched <- rbind(matched, both[, .(spec, n, sim, estimator)])
  }
  expect_equal(data.table::uniqueN(matched), nrow(new))   # every row written has a stored row
})

r2500 <- function(stem, ext = ".csv.gz") stored_path("cv", "our dgp", "r2500", paste0(stem, ext))
versions <- function(g) if (is.list(g$stored)) g$stored else list(g$stored)

test_that("study2.R draws the stored lambda grids: overlap, the three fold counts, both denser specifications", {
  used <- c("overlap", "folds", "denser")
  if (!all(file.exists(r2500(unlist(lapply(study2[used], versions)), "_lamgrid.csv"))))
    skip("lambda grid files not present")
  for (name in used) {
    g <- study2[[name]]
    tmp <- tempfile(fileext = ".csv.gz")
    pilot_grids(g$cells, g$draw, cell_path_settings, tmp)
    new <- data.table::fread(sub("\\.csv\\.gz$", "_lamgrid.csv", tmp))
    by  <- intersect(c("spec", "folds"), names(new))
    for (j in seq_along(versions(g))) {
      version <- names(versions(g))[j]
      rows <- if (is.null(version)) new else new[as.character(get(by)) == version]
      old  <- data.table::fread(r2500(versions(g)[[j]], "_lamgrid.csv"))
      expect_equal(as.data.frame(rows[, .(cell = block, k, lambda)]), as.data.frame(old),
                   tolerance = 1e-12, label = paste(name, version))
    }
  }
})

test_that("study2.R pairs with the stored folds and denser files, one replication of the end cells", {
  used <- c("folds", "denser")
  if (!all(file.exists(r2500(unlist(lapply(study2[used], versions)))))) skip("results files not present")
  tmp <- tempfile("run"); dir.create(tmp)
  withr::local_envvar(c(N_SIM = "1", OUT_DIR = tmp, WORKERS = "1", GRIDS = "^(folds|denser)$",
                        MAX_N = NA, CELLS = "ends"))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(here::here("replication", "runs", "study2.R")),
                    stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0L)
  for (name in used) {
    new <- data.table::fread(file.path(tmp, paste0(name, ".csv.gz")))
    data.table::setnames(new, sub("^exponential", "exp", sub("^quadratic", "quad1", names(new))))
    values <- grep("^(linear|quad1|exp)", names(new), value = TRUE)
    by     <- intersect(c("spec", "folds"), names(new))
    design <- setdiff(names(new), c(values, by, "sim", "estimator", "err"))
    for (version in names(study2[[name]]$stored)) {
      rows <- new[as.character(get(by)) == version]
      old  <- data.table::fread(r2500(study2[[name]]$stored[[version]]))[sim == 1]
      data.table::setnames(old, sub("^exponential", "exp", sub("^quadratic", "quad1", names(old))))
      both <- merge(rows, old, by = c(design, "sim", "estimator"), suffixes = c("", ".stored"))
      expect_equal(nrow(both), nrow(rows), label = paste(name, version, "rows matched"))
      expect_equal(unname(as.matrix(both[, values, with = FALSE])),
                   unname(as.matrix(both[, paste0(values, ".stored"), with = FALSE])),
                   tolerance = 1e-6, label = paste(name, version))
    }
  }
})
