# reruns a few reps of every grid through runs/ on the stored streams and compares w/ stored rows
#   Rscript replication/smoke_test.R   15 reps of the first and last cell of each grid. n up to 50k

library(data.table)

# CELLS=all MAX_N=Inf for every cell incl n = 1 million. SKIP_RUNS=1 OUT_DIR=<folder> recompares
# PARTS=study1 for study 1 only. REPORT=<file> instead of replication/smoke_report.csv
n_sim   <- Sys.getenv("N_SIM", "15")
out_dir <- Sys.getenv("OUT_DIR", file.path(tempdir(), "smoke"))
Sys.setenv(N_SIM = n_sim, MAX_N = Sys.getenv("MAX_N", "50000"), CELLS = Sys.getenv("CELLS", "ends"),
           WORKERS = Sys.getenv("WORKERS", max(1L, parallel::detectCores() - 1L)))

parts <- c(study1 = "study1.R", study2 = "study2.R", wc = "wong_chan.R")
parts <- parts[grepl(Sys.getenv("PARTS", "."), names(parts))]
if (Sys.getenv("SKIP_RUNS") == "")
  for (part in names(parts)) {
    Sys.setenv(OUT_DIR = file.path(out_dir, part))
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(here::here("replication", "runs", parts[[part]])))
    if (status != 0) stop(parts[[part]], " stopped with status ", status)
  }

# the grids w/ the stored files each one reproduces
grids <- lapply(c(study1 = "study1.R", study2 = "study2.R"), function(script)
  local({ source(here::here("replication", "runs", script), local = TRUE); grids }))

stored_file <- function(part, file) {
  where <- switch(part, study1 = list("aipw"), wc = list("cv", "wz_replication"),
                  study2 = list(c("cv/our dgp/r2500", "cv/our dgp/r1000")))
  path <- do.call(file.path, c(list(here::here("results")), where, list(file)))
  path[file.exists(path)][1]
}

#' relative diff so ESS near 1000 and an ATT near -20 are judged on the same scale
rel_diff <- function(a, b) abs(a - b) / pmax(1, abs(b))
#' how far a cell mean moves in MCSEs of the stored run (sd of stored rows / sqrt(stored reps))
mean_shift <- function(a, b, R) {
  d <- abs(mean(a) - mean(b))
  if (is.na(d) || d == 0) 0 else d / (sd(b) / sqrt(R))
}
is_diagnostic <- function(estimator) grepl("^(smd|ess|wmax|lam|nnz)", estimator)
as_key <- function(v) if (is.numeric(v)) format(v, digits = 15, scientific = FALSE, drop0trailing = TRUE, trim = TRUE) else as.character(v)

read_rows <- function(file) {
  x <- fread(file, colClasses = list(character = "err"))
  setnames(x, sub("^exponential", "exp", sub("^quadratic", "quad1", names(x))))
}

#' rerun rows vs a stored file on the cells and estimators that file holds
compare_rows <- function(new, old) {
  values <- grep("^(linear|quad1|exp)", names(new), value = TRUE)
  design <- setdiff(names(new), c(values, "sim", "estimator", "err"))
  for (d in design) { set(new, j = d, value = as_key(new[[d]])); set(old, j = d, value = as_key(old[[d]])) }
  new <- new[unique(old[, design, with = FALSE]), on = design, nomatch = NULL][estimator %in% old$estimator]
  stored_reps <- old[, .(R = uniqueN(sim)), by = design]
  old <- old[unique(new[, c(design, "sim"), with = FALSE]), on = c(design, "sim"), nomatch = NULL]
  timing <- c("time_bal", "time_ps", "time_or")
  new <- new[!estimator %in% timing]
  old <- old[!estimator %in% timing]
  both <- merge(new, old, by = c(design, "sim", "estimator"), suffixes = c("", ".stored"))
  diff <- do.call(cbind, lapply(values, function(v) rel_diff(both[[v]], both[[paste0(v, ".stored")]])))
  na_mismatch <- sum(vapply(values, function(v) sum(xor(is.na(both[[v]]), is.na(both[[paste0(v, ".stored")]]))), 0))
  diag <- is_diagnostic(both$estimator)
  largest <- function(rows) if (any(rows)) max(diff[rows, , drop = FALSE], na.rm = TRUE) else NA_real_
  differing <- both[apply(diff > 1e-10, 1, any, na.rm = TRUE), design, with = FALSE]
  shift <- merge(both, stored_reps, by = design)[, lapply(values, function(v)
             mean_shift(get(v), get(paste0(v, ".stored")), R[1])), by = c(design, "estimator")]
  shift <- unlist(shift[, !c(design, "estimator"), with = FALSE])
  # stored rows not rewritten = fit diagnostics the thesis does not use. counted not compared
  data.table(cells = uniqueN(new[, design, with = FALSE]), replications = uniqueN(new$sim),
             rows_compared = nrow(both), cells_with_diff = uniqueN(differing),
             rows_not_in_stored = nrow(new) - nrow(both),
             largest_diff_estimates = largest(!diag), largest_diff_diagnostics = largest(diag),
             largest_shift_in_mcse = if (all(is.na(shift))) NA_real_ else max(shift, na.rm = TRUE),
             na_mismatch = na_mismatch,
             stored_rows_not_rewritten = uniqueN(old$estimator[!old$estimator %in% new$estimator]))
}

#' rerun lambda grids vs a stored file's. its cell numbers are the blocks
compare_lamgrid <- function(new, old) {
  both <- merge(new[, .(cell = block, k, lambda)], old, by = c("cell", "k"), suffixes = c("", ".stored"))
  data.table(cells = uniqueN(both$cell), replications = NA_integer_, rows_compared = nrow(both),
             cells_with_diff = uniqueN(both$cell[rel_diff(both$lambda, both$lambda.stored) > 1e-10]),
             rows_not_in_stored = nrow(new) - nrow(both),
             largest_diff_estimates = max(rel_diff(both$lambda, both$lambda.stored)),
             largest_diff_diagnostics = NA_real_, na_mismatch = 0L, stored_rows_not_rewritten = 0L)
}

#' one grid file vs every stored file its grid lists. stored files w/ no rerun cell left out
compare_grid <- function(part, new_file) {
  grid <- grids[[part]][[sub("\\.csv\\.gz$", "", basename(new_file))]]
  new  <- read_rows(new_file)
  lam  <- sub("\\.csv\\.gz$", "_lamgrid.csv", new_file)
  lam  <- if (file.exists(lam)) fread(lam)
  by   <- intersect(c("spec", "folds"), names(new))
  stored <- if (is.list(grid$stored)) grid$stored else list(grid$stored)
  rbindlist(lapply(seq_along(stored), function(j) {
    version <- if (is.null(names(stored))) NA_character_ else names(stored)[j]
    keep <- function(x) if (is.na(version)) rep(TRUE, nrow(x)) else as.character(x[[by]]) == version
    rows <- new[keep(new), setdiff(names(new), by), with = FALSE]
    rbindlist(lapply(stored[[j]], function(stem) {
      file <- paste0(stem, ".csv.gz")
      old_file <- stored_file(part, file)
      if (is.na(old_file)) return(data.table(version = version, file = file, stored = FALSE, no_stored_file))
      old <- read_rows(old_file)
      # n = 1 million file is 10 batches of 1000. renumbered as one run of 10000
      if ("batch" %in% names(old)) {
        per <- max(old$sim)
        old[, sim := (batch - 1L) * per + sim][, batch := NULL]
      }
      row <- compare_rows(copy(rows), old)
      out <- if (row$cells > 0) data.table(version = version, file = file, stored = TRUE, row)
      old_lam <- stored_file(part, paste0(stem, "_lamgrid.csv"))
      if (!is.null(lam) && !is.na(old_lam))
        out <- rbind(out, data.table(version = version, file = basename(old_lam), stored = TRUE,
                                     compare_lamgrid(lam[keep(lam)], fread(old_lam))), fill = TRUE)
      out
    }), fill = TRUE)
  }), fill = TRUE)
}

#' one wong chan cell vs its stored .rds
compare_rds <- function(new_file, old_file) {
  new <- readRDS(new_file); old <- readRDS(old_file)
  reps <- seq_len(new$n_rep)
  largest <- function(fields) max(vapply(reps, function(i)
    max(rel_diff(unlist(new$res[[i]][fields]), unlist(old$res[[i]][fields]))), numeric(1)))
  data.table(cells = 1L, replications = length(reps), rows_compared = length(reps), cells_with_diff = NA_integer_,
             rows_not_in_stored = 0L,
             largest_diff_estimates = max(largest(c("est_path", "est_sel", "satt")), rel_diff(new$lam, old$lam)),
             largest_diff_diagnostics = largest(c("lam_end", "lam_sel", "est_mle", "lam_glmnet", "clipped")),
             na_mismatch = 0L, stored_rows_not_rewritten = 0L)
}

no_stored_file <- data.table(cells = NA_integer_, replications = NA_integer_, rows_compared = NA_integer_,
                             cells_with_diff = NA_integer_, rows_not_in_stored = NA_integer_,
                             largest_diff_estimates = NA_real_, largest_diff_diagnostics = NA_real_,
                             largest_shift_in_mcse = NA_real_, na_mismatch = NA_integer_,
                             stored_rows_not_rewritten = NA_integer_)

report <- rbindlist(lapply(names(parts), function(part) {
  files <- list.files(file.path(out_dir, part), "(\\.csv\\.gz|\\.rds)$")
  files <- files[file.size(file.path(out_dir, part, files)) > 0]
  rbindlist(lapply(files, function(f) {
    new <- file.path(out_dir, part, f)
    if (part != "wc") {
      rows <- compare_grid(part, new)
      return(if (nrow(rows)) data.table(part = part, run = f, rows))
    }
    old <- stored_file(part, f)
    data.table(part = part, run = f, version = NA_character_, file = f, stored = !is.na(old),
               if (is.na(old)) no_stored_file else compare_rds(new, old))
  }), fill = TRUE)
}), fill = TRUE)

report[, largest := pmax(largest_diff_estimates, largest_diff_diagnostics, na.rm = TRUE)]
# identical < 1e-10 rel. last digits < 1e-3 as fits near a short path end are rounding sensitive
# means within 0.5 MCSE = no cell mean moves over half its MCSE at the stored rep count
report[, verdict := fcase(!stored, "no stored file here", na_mismatch > 0, "check",
                          largest < 1e-10, "identical", largest < 1e-3, "last digits",
                          largest_shift_in_mcse < 0.5, "means within 0.5 MCSE", default = "differs")]
fwrite(report, Sys.getenv("REPORT", here::here("replication", "smoke_report.csv")))
print(report[, .(run, version, file, cells, replications, rows_compared, cells_with_diff, largest = signif(largest, 2),
                 shift_in_mcse = signif(largest_shift_in_mcse, 2), verdict)])
