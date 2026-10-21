# recomputes the study 1 section summaries and coverage_empse.rds from results/aipw w/ R/measures.R
#   Rscript replication/verify_summaries.R   about 7 minutes. needs results/ (7 GB)

source(here::here("replication", "R", "measures.R"))
stored_outcome <- c(linear = "linear", quadratic = "quad1", exponential = "exp")
summaries <- list.files(here::here("output", "thesis_summaries"), "\\.rds$", full.names = TRUE)

#' largest diff relative to the stored number where that is over 1
#' relative as big sums of squares lose last digits w/o extended precision (apple silicon)
largest <- function(a, b) {
  if (length(a) && length(b)) max(abs(a - b) / pmax(1, abs(b)), na.rm = TRUE) else NA_real_
}
cat("R", R.version$major, R.version$minor, R.version$platform, "long double bytes",
    .Machine$sizeof.longdouble, "\n\n")

report <- list()
not_here <- character()
for (f in summaries) {
  s <- readRDS(f)
  id <- sub("\\.rds$", "", basename(f))
  files <- s$info$files
  # only study 1 sections list files w/ a spec. coverage_empse is checked below
  if (is.null(files) || !all(c("spec", "file", "arms") %in% names(files)) || is.null(s$cells)) {
    if (id != "coverage_empse") not_here <- c(not_here, id)
    next
  }
  paths <- here::here("results", "aipw", files$file)
  if (!all(file.exists(paths))) { report[[id]] <- "skipped: results file absent"; next }

  results <- Map(function(path, arms, spec) {
    r <- read_results(path, trimws(strsplit(arms, ",")[[1]]))
    mutate(r, spec = spec)
  }, paths, files$arms, files$spec) |> bind_rows()
  # add-on file repeats oracle_hajek so keep it once
  # a summary using part of a grid keeps the levels its stored table holds
  results <- distinct(results, spec, level, sim, estimator, outcome, .keep_all = TRUE) |>
    semi_join(distinct(as_tibble(s$perf), spec, level), by = c("spec", "level"))

  perf <- performance(results) |> mutate(outcome = stored_outcome[outcome])
  j <- inner_join(perf, as_tibble(s$perf), by = c("spec", "level", "estimator", "outcome"),
                  suffix = c("", ".stored"))
  diffs <- c(perf_rows_matched = nrow(j) - nrow(s$perf),
             bias = largest(j$bias, j$bias.stored), empse = largest(j$empse, j$empse.stored),
             rmse = largest(j$rmse, j$rmse.stored), mcse_rmse = largest(j$mcse_rmse, j$mcse_rmse.stored))
  for (i in seq_len(nrow(s$info$pairs))) {
    A <- s$info$pairs$A[i]; B <- s$info$pairs$B[i]
    if (!all(c(A, B) %in% results$estimator)) next
    pd <- pair_draws(results, A, B)
    cells <- pair_cells(pd$draws) |> mutate(outcome = stored_outcome[outcome])
    sc <- as_tibble(s$cells) |> filter(A == !!A, B == !!B, kept)
    jc <- inner_join(cells, sc, by = c("spec", "level", "outcome"), suffix = c("", ".stored"))
    avg <- relative_rmse(results, A, B) |>
      mutate(outcome = ifelse(outcome == "pooled", "pooled", stored_outcome[outcome]))
    sa <- as_tibble(s$avg) |> filter(A == !!A, B == !!B, theta == "pair", set == "all", is.na(factor))
    ja <- inner_join(avg, sa, by = c("block", "outcome"), suffix = c("", ".stored"))
    diffs <- c(diffs, setNames(c(nrow(jc) - nrow(sc), largest(jc$R_c, jc$R_c.stored),
                                 largest(jc$log_R_c_mcse, jc$log_R_c_mcse.stored),
                                 largest(ja$R_bar, ja$R_bar.stored), largest(ja$R_bar_mcse, ja$R_bar_mcse.stored)),
                               paste(A, "vs", B, c("cells_matched", "R_c", "R_c_mcse", "R_bar", "R_bar_mcse"))))
  }
  report[[id]] <- diffs
  cat(sprintf("%-30s largest difference %.1e\n", id, max(diffs[!grepl("matched", names(diffs))], na.rm = TRUE)))
}

# coverage at empirical SE (Fig 6 and App Figs 2 and 5) and the median MLE weight diagnostics
# built as the qmd chunk cov-data builds them
cov_files <- list(
  n       = list(correct = "aipw_correct_n400k_dgp2.csv.gz",
                 misspec = c("aipw_misspec_n_dgp2.csv.gz", "aipw_misspec_n_dgp2n1mill.csv.gz")),
  overlap = list(correct = c("aipw_correct_overlap7_dgp2.csv.gz", "aipw_correct_overlap8_dgp2.csv.gz"),
                 misspec = c("aipw_misspec_overlap7_dgp2.csv.gz", "aipw_misspec_overlap8_dgp2.csv.gz")),
  signal  = list(correct = "aipw_correct_strength_dgp2.csv.gz",
                 misspec = "aipw_misspec_strength_dgp2.csv.gz"))
cov_arms  <- c("balnet0", "abw_balnet0", "glmnetcv_hajek", "aipw_glmnetcv", "oracle_hajek")
diag_arms <- c("wmax_glm1", "ess_glm1")
stored_cov <- readRDS(here::here("output", "thesis_summaries", "coverage_empse.rds"))
for (grid in names(cov_files)) for (spec in names(cov_files[[grid]])) {
  id <- paste("coverage_empse", grid, spec)
  paths <- here::here("results", "aipw", cov_files[[grid]][[spec]])
  if (!all(file.exists(paths))) { report[[id]] <- "skipped: results file absent"; next }
  arms <- if (grid == "overlap") c(cov_arms, diag_arms) else cov_arms
  results <- bind_rows(lapply(paths, read_results, estimators = arms))
  keys <- c(intersect(c("n", "overlap", "strength"), names(stored_cov[[grid]])), "estimator", "outcome")
  perf <- performance(filter(results, estimator %in% cov_arms)) |> mutate(outcome = stored_outcome[outcome])
  sc <- as_tibble(stored_cov[[grid]]) |> filter(spec == !!spec)
  j <- inner_join(perf, sc, by = keys, suffix = c("", ".stored"))
  diffs <- c(rows_matched = nrow(j) - nrow(sc), m = largest(j$m, j$m.stored),
             coverage = largest(j$coverage, j$coverage.stored), mcse = largest(j$mcse_coverage, j$mcse),
             bias = largest(j$bias, j$bias.stored), empse = largest(j$empse, j$empse.stored))
  if (grid == "overlap") {
    wd <- results |>
      filter(estimator %in% diag_arms, outcome == "linear", is.finite(error)) |>
      summarise(value = median(error), .by = c(n, overlap, estimator))
    sw <- as_tibble(stored_cov$diag) |> filter(spec == !!spec)
    jw <- inner_join(wd, sw, by = c("n", "overlap", "estimator"), suffix = c("", ".stored"))
    diffs <- c(diffs, diag_rows_matched = nrow(jw) - nrow(sw), diag_median = largest(jw$value, jw$value.stored))
  }
  report[[id]] <- diffs
  cat(sprintf("%-30s largest difference %.1e\n", id, max(diffs[!grepl("matched", names(diffs))], na.rm = TRUE)))
}

cat("\nDetail\n")
for (id in names(report)) {
  cat("\n==", id, "\n")
  if (is.character(report[[id]])) cat("  ", report[[id]], "\n") else print(signif(report[[id]], 3))
}
cat("\nNot recomputed here, checked by handoff/scripts/build_summaries.R (verification.md):\n",
    paste(not_here, collapse = " "), "\n")

# stop if any diff is over 1e-10 or a stored row is not recomputed
failed <- names(report)[vapply(report, function(d) is.numeric(d) &&
  (any(d[grepl("matched", names(d))] != 0) || any(d[!grepl("matched", names(d))] > 1e-10, na.rm = TRUE)),
  logical(1))]
if (length(failed)) stop("differences above 1e-10 or unmatched rows in: ", paste(failed, collapse = ", "))
