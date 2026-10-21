# stored data for the appendix full tables (perf diagnostics approx bias vs OLS and 800k ipw run)
# build_summaries.R saves output/thesis_summaries/<id>.rds. coverage in sections.R

library(dplyr)
library(tidyr)
library(purrr)

#' x or y if x is NULL. base R only has it from 4.4
`%||%` <- function(x, y) if (is.null(x)) y else x

# read

#' some arms of one results file. one row per draw and arm w/ outcome models as cols
#' grid_key pasted as eg 4e+05 so the MAE mcse bootstrap hits the same cells
read_wide <- function(file, arms, dir = "aipw") {
  d <- data.table::fread(here::here("results", dir, file),
                         colClasses = list(character = "err"))
  d <- d[d$estimator %in% arms, ]
  d$err <- NULL
  oc <- intersect(outcome_names, names(d))
  for (o in oc) d[[o]] <- suppressWarnings(as.numeric(d[[o]]))
  design <- setdiff(names(d), c("sim", "estimator", oc))
  as_readr_text <- \(v) if (is.numeric(v)) as.character(as.numeric(v)) else
    as.character(v)
  d$level <- do.call(paste, c(lapply(d[, ..design], key_text), sep = "_"))
  d$grid_key <- do.call(paste, c(lapply(d[, ..design], as_readr_text),
                                 sep = "_"))
  tibble::as_tibble(d) |>
    relocate(level, grid_key)
}

#' int design cols to double like readr gave them so values format as before
as_doubles <- function(design) {
  mutate(design, across(where(is.integer), as.numeric))
}

#' every file of one spec of an appendix entry wide (read_wide()) plus derived design cols
read_spec_wide <- function(entry, spec, arms = NULL) {
  wide <- imap(entry$files[[spec]], \(listed, file) {
    read_wide(file, arms %||% listed, entry$dir %||% "aipw") |>
      mutate(file = file, .before = 1)
  }) |>
    list_rbind()
  if (!is.null(entry$derive)) wide <- entry$derive(wide)
  wide
}

# appendix table data

#' full perf of every arm and cell of one spec. cells in grid_key order after set.seed(seed)
#' caller's rng stream restored after
full_performance <- function(wide, B, seed, miss_max = 0.01) {
  oc <- intersect(outcome_names, names(wide))
  design <- wide |>
    select(-file, -grid_key, -sim, -estimator, -all_of(oc)) |>
    distinct() |>
    as_doubles()
  keys <- distinct(wide, grid_key, level)
  stream <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  on.exit(if (is.null(stream)) {
    rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", stream, envir = globalenv())
  })
  set.seed(seed)
  wide |>
    select(grid_key, estimator, all_of(oc)) |>
    pivot_longer(all_of(oc), names_to = "outcome", values_to = "error") |>
    group_by(grid_key, estimator, outcome) |>
    summarise(arm_measures_full(error, B), .groups = "drop") |>
    mutate(dropped = missing / reps > miss_max) |>
    left_join(keys, by = "grid_key") |>
    select(-grid_key) |>
    relocate(level) |>
    left_join(design, by = "level")
}

#' fit diagnostics per cell over reps. trunc counts paths ending above trunc_at (balance tol)
#' only nnz_or* differ by outcome. rest is logged once per fit so kept from the first col
diag_summaries <- function(wide, trunc_at = 1.01e-4) {
  oc <- intersect(outcome_names, names(wide))
  design <- wide |>
    select(-file, -grid_key, -sim, -estimator, -all_of(oc)) |>
    distinct() |>
    as_doubles()
  repeated <- as.matrix(select(filter(wide, !startsWith(estimator, "nnz_or")),
                               all_of(oc)))
  agree <- ncol(repeated) < 2 || isTRUE(all(
    repeated == repeated[, 1] | (is.na(repeated) & is.na(repeated[, 1]))))
  diag <- wide |>
    select(level, estimator, all_of(oc)) |>
    pivot_longer(all_of(oc), names_to = "outcome", values_to = "v") |>
    filter(startsWith(estimator, "nnz_or") | outcome == oc[1]) |>
    mutate(outcome = if_else(startsWith(estimator, "nnz_or"), outcome, NA)) |>
    summarise(reps = n(), missing = sum(!is.finite(v)),
              centre = median(v, na.rm = TRUE),
              lo = quantile(v, 0.25, na.rm = TRUE, names = FALSE),
              hi = quantile(v, 0.75, na.rm = TRUE, names = FALSE),
              sd = sd(v, na.rm = TRUE), min = min(v, na.rm = TRUE),
              max = max(v, na.rm = TRUE),
              .by = c(level, estimator, outcome)) |>
    mutate(stat = "median") |>
    left_join(design, by = "level")
  trunc <- wide |>
    filter(estimator %in% c("lam_end1", "lam_end0")) |>
    summarise(k = sum(.data[[oc[1]]] > trunc_at, na.rm = TRUE),
              reps = n(), .by = c(level, estimator))
  list(diag = diag, trunc = trunc, outcomes_agree = agree)
}

#' os R version and run date from a results file's _session.txt (written by R/simulate.R)
#' for the fit times note. NULL if there is none
session_text <- function(file, dir = "aipw") {
  f <- here::here("results", dir, sub("\\.csv\\.gz$", "_session.txt", file))
  if (!file.exists(f)) return(NULL)
  s <- readLines(f, warn = FALSE)
  pick <- \(pattern) sub(pattern, "", grep(pattern, s, value = TRUE)[1])
  r <- sub(" \\(.*", "", pick("^R version "))
  os <- pick("^Running under: ")
  d <- substr(pick("^date: "), 1, 10)
  paste0(os, ", R ", r, if (!is.na(d)) paste0(", run ", d))
}

# one appendix entry

#' everything one appendix entry's tables show from its results files. missing parts NULL
build_appendix <- function(entry, id, miss_max = 0.01) {
  dir <- entry$dir %||% "aipw"
  specs <- names(entry$files)
  checks <- list()
  add_check <- \(check, ok, detail) {
    checks[[length(checks) + 1]] <<- check_row(check, ok, detail)
  }
  out <- list()

  if (!is.null(entry$full)) {
    out$full <- map(set_names(specs), \(spec) {
      wide <- read_spec_wide(entry, spec)
      wide <- filter(wide, !estimator %in% entry$full$drop)
      # an arm listed for two files would count its draws twice
      dup <- anyDuplicated(select(wide, level, sim, estimator))
      add_check(paste0("full (", spec, "): each arm read from one file"),
                dup == 0, if (dup == 0) "no duplicated (level, sim, arm)" else
                  "duplicated draws")
      perf <- full_performance(wide, entry$full$B, entry$full$seed, miss_max)
      add_check(paste0("full (", spec, "): replications per arm and cell = ",
                       entry$reps), all(perf$reps == entry$reps),
                sprintf("min %d, max %d", min(perf$reps), max(perf$reps)))
      mutate(perf, spec = spec, .before = 1)
    })
  }

  if (isTRUE(entry$diag)) {
    d <- map(set_names(specs), \(spec) {
      s <- diag_summaries(read_spec_wide(entry, spec, arms_diag))
      add_check(paste0("diagnostics (", spec, "): logged once per fit, ",
                       "the same in every outcome column"), s$outcomes_agree,
                "every diagnostic but nnz_or* agrees across outcome columns")
      add_check(paste0("diagnostics (", spec, "): replications per cell = ",
                       entry$reps), all(s$diag$reps == entry$reps),
                sprintf("min %d, max %d", min(s$diag$reps), max(s$diag$reps)))
      s
    })
    out$diag <- map(d, "diag")
    out$trunc <- map(d, "trunc")
    out$session <- map(set_names(specs), \(spec) {
      session_text(names(entry$files[[spec]])[1], dir)
    })
  }

  if (!is.null(entry$approx)) {
    A <- entry$approx$A
    B <- entry$approx$B
    out$approx <- map(set_names(specs), \(spec) {
      files <- entry$files[[spec]]
      # check arm read from every file only to test the pairing below
      est <- imap(files, \(listed, file) {
        read_results(file, union(intersect(listed, c(A, B)), entry$check_arm),
                     dir) |>
          mutate(file = file)
      }) |>
        list_rbind() |>
        mutate(spec = spec, .before = 1)
      # A and B from different files pair by level and sim so the files must hold the same draws
      # check arm is in every file and must agree rep by rep
      if (!is.null(entry$check_arm) && length(files) > 1) {
        first <- filter(est, file == names(files)[1],
                        estimator == entry$check_arm)
        for (this_file in names(files)[-1]) {
          cmp <- compare_rows(first, filter(est, file == this_file,
                                            estimator == entry$check_arm))
          add_check(paste0("add-on pairing (", spec, "): ", entry$check_arm,
                           " in ", this_file, " = ", names(files)[1]),
                    cmp$matched == cmp$n_x && cmp$matched == cmp$n_y &&
                      cmp$max_diff <= 1e-10,
                    sprintf("%d of %d rows matched, max |diff| %.2g",
                            cmp$matched, cmp$n_x, cmp$max_diff))
        }
      }
      # each arm only from the file it is listed for
      listed <- imap(files, \(arms, file) {
        tibble(file = file, estimator = intersect(arms, c(A, B)))
      }) |>
        list_rbind()
      est <- est |>
        semi_join(listed, by = c("file", "estimator")) |>
        select(-file)
      if (!is.null(entry$derive)) est <- entry$derive(est)
      design <- est |>
        select(-spec, -sim, -estimator, -outcome, -error) |>
        distinct() |>
        as_doubles()
      pd <- pair_draws(est, A, B, miss_max)
      n_cells <- nrow(pd$cells)
      kept <- sum(pd$cells$kept)
      add_check(paste0("approximation bias (", spec, "): ", A, " and ", B,
                       " on the same draws, each cell once"),
                anyDuplicated(select(pd$draws, level, sim, outcome)) == 0 &&
                  all(count(pd$draws, level, outcome)$n <= entry$reps),
                sprintf("%d of %d cells kept", kept, n_cells))
      list(cells = approx_cell_measures(pd$draws) |>
             left_join(design, by = "level") |>
             mutate(A = A, B = B),
           avg = mutate(approx_average(pd$draws), A = A, B = B),
           n_levels = n_distinct(est$level), n_cells = n_cells,
           kept = kept)
    })
  }

  if (isTRUE(entry$simple)) {
    est <- imap(entry$files$correct, \(listed, file) {
      read_results(file, listed, dir)
    }) |>
      list_rbind()
    design <- est |>
      select(-sim, -estimator, -outcome, -error) |>
      distinct() |>
      as_doubles()
    dup <- anyDuplicated(select(est, level, sim, estimator, outcome))
    add_check("each draw read once", dup == 0,
              if (dup == 0) "no duplicated (level, sim, arm, outcome)" else
                "duplicated draws")
    out$simple <- est |>
      summarise(arm_measures(error), .by = c(level, estimator, outcome)) |>
      left_join(design, by = "level")
    add_check(paste("replications per arm and cell =", entry$reps),
              all(out$simple$reps == entry$reps) &&
                all(out$simple$missing == 0),
              sprintf("%d rows; min %d, max %d replications, %d missing",
                      nrow(out$simple), min(out$simple$reps),
                      max(out$simple$reps), sum(out$simple$missing)))
  }

  design_cols <- function() {
    first <- out$full[[1]] %||% out$diag[[1]] %||% out$simple
    setdiff(names(first), c(
      "spec", "level", "estimator", "outcome", "reps", "missing", "m",
      "n_extreme", "bias", "empse", "rmse", "mcse_bias", "mcse_empse",
      "mcse_rmse", "share", "mcse_share", "kurt", "mcse_empse_k", "mae",
      "mcse_mae", "dropped", "centre", "lo", "hi", "sd", "min", "max",
      "stat"))
  }
  c(list(info = list(id = id, title = entry$title, kind = "appendix",
                     dir = dir, specs = specs,
                     files = section_files(entry$files, dir),
                     design = design_cols(), display = entry$display,
                     reps = entry$reps, miss_max = miss_max, tau = 0,
                     B = entry$full$B, seed = entry$full$seed,
                     built = Sys.time(), r_version = R.version.string,
                     packages = map_chr(c("dplyr", "tidyr", "data.table"),
                                        \(p) as.character(packageVersion(p)))),
         checks = list_rbind(checks)),
    out)
}
