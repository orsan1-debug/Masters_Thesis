# results files to one stored summary per section. reads + checks pairing + computes what it shows
# the measures live in rbar.R. this only arranges the data and picks what is averaged over what

library(dplyr)
library(tidyr)
library(purrr)

outcome_names <- c("linear", "quad1", "exp")

# estimator rows that are logged diagnostics (path ends/ESS/timings) not estimates
diag_pattern <- paste0("^(lam_|nnz_|smd|ess|wmax|time_|cvloss_|logcvloss_|",
                       "mse_|msre_|mu_|trunc05|prev|emin|emax|nout)")

# read

#' design value as text whatever type fread gave so 100000 and 1e+05 both give "100000"
key_text <- function(v) {
  if (is.numeric(v)) {
    format(as.numeric(v), scientific = FALSE, drop0trailing = TRUE,
           trim = TRUE, digits = 15)
  } else {
    as.character(v)
  }
}

#' rows of some arms from one results file. design cols pasted into `level` so steps group on it
read_results <- function(file, arms, dir = "aipw") {
  path <- here::here("results", dir, file)
  d <- data.table::fread(path, colClasses = list(character = "err"))
  d <- d[d$estimator %in% arms, ]
  d$err <- NULL                           # the runner's error messages
# n = 1M run is one cell in 10 batches (own checkpoints in simulate_grid())
# number sims across batches so the rep check and MCSEs count them together
  if ("batch" %in% names(d)) {
    d$sim <- (d$batch - 1L) * max(d$sim) + d$sim
    d$batch <- NULL
  }
  oc <- intersect(outcome_names, names(d))
  design <- setdiff(names(d), c("sim", "estimator", oc))
  d$level <- do.call(paste, c(lapply(d[, ..design], key_text), sep = "_"))
  long <- data.table::melt(d, measure.vars = oc, variable.name = "outcome",
                           value.name = "error", variable.factor = FALSE)
# tau = 0 everywhere so the estimate is the error
# cols read as text (failed reps) become numeric w/ NA where it failed
  long$error <- suppressWarnings(as.numeric(long$error))
  tibble::as_tibble(long)
}

#' every file of a section tagged w/ spec and file
#' check arm read from all files so check_section() can compare them. kept from the first only
read_section <- function(section) {
  imap(section$files, \(spec_files, spec) {
    imap(spec_files, \(arms, file) {
      d <- read_results(file, union(arms, section$check_arm))
      if (!is.null(section$keep)) d <- d[section$keep(d), ]
      mutate(d, spec = spec, file = file, .before = 1)
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}

# checks. one row each. build_summaries.R stores nothing if one fails
# as every average assumes both arms of a pair used the same draws

#' one row of the checks table
check_row <- function(check, ok, detail) {
  tibble(check = check, ok = isTRUE(ok), detail = detail)
}

#' compare one arm in two tables draw by draw. rows matched + largest abs diff
compare_rows <- function(x, y) {
  key <- c("level", "sim", "outcome")
  j <- inner_join(select(x, all_of(key), ex = error),
                  select(y, all_of(key), ey = error), by = key)
  list(n_x = nrow(x), n_y = nrow(y), matched = nrow(j),
       max_diff = suppressWarnings(max(abs(j$ex - j$ey), na.rm = TRUE)))
}

#' check the files pair as designed then keep each arm only from its listed file
check_section <- function(res, section) {
  checks <- list()
  tol <- 1e-10

# add-on files add arms on the same draws as the main file
# only true if the check arm in both agrees rep by rep
  for (this_spec in names(section$files)) {
    files <- names(section$files[[this_spec]])
    if (length(files) < 2 || is.null(section$check_arm)) next
    first <- filter(res, spec == this_spec, file == files[1],
                    estimator == section$check_arm)
    for (this_file in files[-1]) {
      other <- filter(res, spec == this_spec, file == this_file,
                      estimator == section$check_arm)
      cmp <- compare_rows(first, other)
      checks[[length(checks) + 1]] <- check_row(
        paste0("add-on pairing (", this_spec, "): ", section$check_arm,
               " in ", this_file, " = ", files[1]),
        cmp$matched == cmp$n_x && cmp$matched == cmp$n_y &&
          cmp$max_diff <= tol,
        sprintf("%d of %d rows matched, max |diff| %.2g", cmp$matched,
                cmp$n_x, cmp$max_diff))
    }
  }
# done w/ the check arm copies. keep each arm only from the file it is listed for
  listed <- imap(section$files, \(spec_files, this_spec) {
    imap(spec_files, \(arms, this_file) {
      tibble(spec = this_spec, file = this_file, estimator = arms)
    }) |>
      list_rbind()
  }) |>
    list_rbind()
  res <- semi_join(res, listed, by = c("spec", "file", "estimator"))

# no arm read from two files or a draw counts twice
  n_unique <- nrow(distinct(res, spec, level, sim, outcome, estimator))
  checks[[length(checks) + 1]] <- check_row(
    "each arm read from one file only", n_unique == nrow(res),
    if (n_unique == nrow(res)) {
      "no duplicated (spec, level, sim, outcome, estimator)"
    } else {
      "duplicated rows: an arm is listed for two files"
    })

# stated reps for every arm in every cell (spec level outcome)
  n_rep <- count(res, spec, level, outcome, estimator)$n
  checks[[length(checks) + 1]] <- check_row(
    paste("replications per arm and cell =", section$reps),
    is.null(section$reps) || all(n_rep == section$reps),
    sprintf("min %d, max %d", min(n_rep), max(n_rep)))

# design cols constant within a level
  design <- setdiff(names(res), c("spec", "file", "level", "sim",
                                  "estimator", "outcome", "error"))
  n_design <- nrow(distinct(res[c("spec", "level", design)]))
  n_level <- nrow(distinct(res[c("spec", "level")]))
  checks[[length(checks) + 1]] <- check_row(
    "one design row per level", n_design == n_level,
    sprintf("%d levels, %d design rows", n_level, n_design))

# "both" block and its MCSE assume both specs ran on the same draws
# so same levels and draws and the oracles (no ps spec) must agree rep by rep
  specs <- names(section$files)
  if (length(specs) > 1) {
# local names must not clash w/ a design col (eg s)
    draws <- map(specs, \(this_spec) {
      distinct(filter(res, spec == this_spec), level, sim)
    })
    same <- nrow(dplyr::setdiff(draws[[1]], draws[[2]])) == 0 &&
      nrow(dplyr::setdiff(draws[[2]], draws[[1]])) == 0
    checks[[length(checks) + 1]] <- check_row(
      "pooled specifications have the same levels and draws", same,
      sprintf("%d and %d draws", nrow(draws[[1]]), nrow(draws[[2]])))
    for (this_arm in section$spec_arms) {
      cmp <- compare_rows(
        filter(res, spec == specs[1], estimator == this_arm),
        filter(res, spec == specs[2], estimator == this_arm))
      checks[[length(checks) + 1]] <- check_row(
        paste0("shared draws across specifications: ", this_arm),
        cmp$matched == cmp$n_x && cmp$matched == cmp$n_y &&
          cmp$max_diff <= tol,
        sprintf("%d of %d rows matched, max |diff| %.2g", cmp$matched,
                cmp$n_x, cmp$max_diff))
    }
  }

# an Inf estimate would silently count as missing (is.finite)
  n_inf <- sum(is.infinite(res$error))
  checks[[length(checks) + 1]] <- check_row(
    "no infinite estimates", n_inf == 0, sprintf("%d infinite", n_inf))

  list(res = select(res, -file), checks = list_rbind(checks))
}

# averages for every block and set

#' design value as text for tables eg 0.5 or 1000
fmt_value <- function(v) {
  if (is.numeric(v)) {
    format(v, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  } else {
    as.character(v)
  }
}

#' spec blocks to average over. each spec alone and "both" when there are two
spec_blocks <- function(specs) {
  blocks <- set_names(as.list(specs), specs)
  if (length(specs) > 1) blocks$both <- specs
  blocks
}

#' Rbar and sbar of one pair over a group of cells. pooled then each outcome
average_group <- function(by_level, theta_all) {
  oc <- intersect(outcome_names, theta_all$outcome)
  map(c("pooled", oc), \(o) {
    group <- if (o == "pooled") theta_all else filter(theta_all, outcome == o)
    used <- filter(group, use)
    average_over_cells(by_level, used) |>
      mutate(outcome = o, n_dropped = nrow(group) - nrow(used), .before = 1)
  }) |>
    list_rbind()
}

#' every average of one pair. each block and set (or level of its by factors) and both thetas
#' theta pair = cells this pair keeps. common = cells every compared arm keeps
average_sets <- function(by_level, cells, blocks, sets) {
  imap(blocks, \(specs, block) {
# a pair only enters a block if both arms ran under all its specs
# an arm never run there is absent not missing
    if (!all(specs %in% cells$spec)) return(NULL)
    imap(sets, \(set, name) {
      theta <- filter(cells, spec %in% specs)
      if (!is.null(set$keep)) theta <- theta[set$keep(theta), ]
      groups <- if (is.null(set$by)) {
        list(list(factor = NA_character_, value = NA_character_,
                  rows = rep(TRUE, nrow(theta))))
      } else {
        list_flatten(map(set$by, \(f) map(sort(unique(theta[[f]])), \(v) {
          list(factor = f, value = fmt_value(v), rows = theta[[f]] == v)
        })))
      }
      map(groups, \(g) {
        map(c("pair", "common"), \(th) {
          t_all <- mutate(theta[g$rows, ],
                          use = kept & (th == "pair" | coalesce(common, FALSE)))
          average_group(by_level, t_all) |>
            mutate(theta = th, factor = g$factor, value = g$value,
                   .before = 1)
        }) |>
          list_rbind()
      }) |>
        list_rbind()
    }) |>
      list_rbind(names_to = "set")
  }) |>
    compact() |>
    list_rbind(names_to = "block")
}

# one section

#' everything one thesis section shows computed from its results files
build_section <- function(section, id, miss_max = 0.01) {
  read <- check_section(read_section(section), section)
  res <- read$res
  if (!is.null(section$derive)) res <- section$derive(res)
  extra <- map(section$extra, \(f) f(res))

# estimates only. logged diagnostics are not arms
  diag_rows <- grep(diag_pattern, unique(res$estimator), value = TRUE)
  est <- filter(res, !estimator %in% diag_rows)
  rm(res)
  design <- est |>
    select(-sim, -estimator, -outcome, -error) |>
    distinct()
  perf <- cell_performance(est, miss_max) |>
    left_join(design, by = c("spec", "level"))

# common Theta = cells where every compared arm under that spec has <= miss_max missing
  pairs <- filter(section$pairs, A %in% est$estimator, B %in% est$estimator)
  compared <- unique(c(pairs$A, pairs$B))
  common <- perf |>
    filter(estimator %in% compared) |>
    summarise(common = !any(dropped), .by = c(spec, level, outcome))

  blocks <- spec_blocks(names(section$files))
  sets <- c(list(all = list()), section$sets)
  by_pair <- pmap(pairs, \(A, B) {
# only specs where both arms were run
    run <- est |>
      filter(estimator %in% c(A, B)) |>
      distinct(spec, estimator) |>
      count(spec) |>
      filter(n == 2)
    est_ab <- semi_join(est, run, by = "spec")
    pd <- pair_draws(est_ab, A, B, miss_max)
    cells <- pd$cells |>
      left_join(pair_cell_measures(pd$draws),
                by = c("spec", "level", "outcome")) |>
      left_join(common, by = c("spec", "level", "outcome")) |>
      left_join(design, by = c("spec", "level"))
    list(cells = mutate(cells, A = A, B = B, .before = 1),
         avg = average_sets(level_matrices(pd$draws), cells, blocks,
                            sets) |>
           mutate(A = A, B = B, .before = 1))
  })

  list(
    info = list(id = id, title = section$title, x = section$x,
                by = section$by, specs = names(section$files),
                files = section_files(section$files), pairs = pairs,
                sets = names(sets), reps = section$reps,
                miss_max = miss_max, tau = 0, built = Sys.time(),
                r_version = R.version.string,
                packages = map_chr(c("dplyr", "tidyr", "data.table"),
                                   \(p) as.character(packageVersion(p)))),
    checks = read$checks,
    perf = perf,
    cells = list_rbind(map(by_pair, "cells")),
    avg = list_rbind(map(by_pair, "avg")),
    extra = extra
  )
}

#' files read w/ size mtime and md5 so read_summary() (thesis.R) can tell if one changed
section_files <- function(files, dir = "aipw") {
  imap(files, \(spec_files, spec) {
    tibble(spec = spec, file = names(spec_files),
           arms = map_chr(spec_files, \(a) paste(a, collapse = ", ")))
  }) |>
    list_rbind() |>
    mutate(path = here::here("results", dir, file),
           size = file.size(path),
           mtime = file.mtime(path),
           md5 = unname(tools::md5sum(path))) |>
    select(-path)
}
