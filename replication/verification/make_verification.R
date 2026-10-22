# writes replication/VERIFICATION.md from run_checks.R logs smoke_report.csv and compare_render.R
#   Rscript replication/verification/make_verification.R   after run_checks.R and compare_render.R

rp <- function(...) here::here("replication", ...)
lines <- function(path) if (file.exists(path)) readLines(path, warn = FALSE) else character()
items <- function(x) trimws(unlist(strsplit(x, ";")))
count <- function(v, levels) {
  n <- table(factor(v, levels))
  paste(n[n > 0], names(n)[n > 0], collapse = ", ")
}
verdicts <- c("identical", "last digits", "means within 0.5 MCSE", "differs")

displays <- read.csv(rp("verification", "displays.csv"), check.names = FALSE)
session  <- lines(rp("sessionInfo.txt"))
steps    <- lines(rp("logs", "run_checks.txt"))
ran_at   <- format(file.mtime(rp("logs", "run_checks.txt")), "%d %b %Y %H:%M")

# code vs the code that ran
test_tally <- function(file) {
  log <- lines(rp("logs", file))
  t <- regmatches(log, regexec("FAIL (\\d+) \\| WARN (\\d+) \\| SKIP (\\d+) \\| PASS (\\d+)", log))
  t <- do.call(rbind, lapply(t[lengths(t) == 5], \(m) as.integer(m[-1])))
  c(files = length(grep("^== test_.*\\.R ==$", log)), pass = sum(t[, 4]), fail = sum(t[, 1]), skip = sum(t[, 3]))
}
tally <- test_tally("tests.txt")
# and in a copy w/o results/ if that run is logged
bare <- if (file.exists(rp("logs", "tests_without_results.txt"))) test_tally("tests_without_results.txt")

# stored summaries vs stored results
summ <- lines(rp("logs", "summaries.txt"))
hit  <- regmatches(summ, regexec("^(\\S.*?)\\s+largest difference (\\S+)$", summ))
hit  <- do.call(rbind, hit[lengths(hit) == 3])
recomputed <- setNames(suppressWarnings(as.numeric(hit[, 3])), trimws(hit[, 2]))
second <- lines(here::here("output", "thesis_summaries", "verification.md"))
second <- regmatches(second, regexec("^## (\\S+) \\((passed|failed)\\)", second))
second <- do.call(rbind, second[lengths(second) == 3])
second <- setNames(second[, 3], second[, 2])

summary_check <- function(ids) {
  if (grepl("^design", ids)) return(sub("^design: ", "", ids))
  ids <- items(ids)
  key <- lapply(ids, \(id) if (grepl("^coverage_empse:", id))
    paste("coverage_empse", sub(".*:", "", id), c("correct", "misspec")) else id)
  here_ <- vapply(key, \(k) all(k %in% names(recomputed)), logical(1))
  out <- character()
  if (any(here_)) {
    d <- vapply(key[here_], \(k) max(recomputed[k]), numeric(1))
    out <- c(out, paste0(ifelse(is.finite(d) & d <= 1e-10, "", "DIFFERS "), ids[here_], " ",
                         format(signif(d, 2)), collapse = ", "))
  }
  two <- ids[!here_ & ids %in% names(second)]
  if (length(two)) out <- c(out, sprintf("%s %s (second computation)", paste(two, collapse = ", "),
                                         if (all(second[two] == "passed")) "passed" else "FAILED"))
  none <- ids[!here_ & !ids %in% names(second)]
  if (length(none)) out <- c(out, paste("NOT CHECKED:", paste(none, collapse = ", ")))
  paste(out, collapse = ". ")
}

# stored results vs a rerun of the code
smoke <- read.csv(rp("smoke_report.csv"))
rerun_check <- function(results) {
  if (!nzchar(results)) return("")
  stems <- sub("\\{.*", "", items(results))
  hits  <- lapply(stems, \(s) startsWith(smoke$file, s) &
                    grepl("^(\\.|_lamgrid|\\d)", substring(smoke$file, nchar(s) + 1)))
  out <- count(smoke$verdict[Reduce(`|`, hits)], verdicts)
  missing <- stems[!vapply(hits, any, logical(1))]
  if (length(missing)) out <- paste0(out, if (nzchar(out)) ". ", "Not rerun: ", paste(missing, collapse = ", "))
  out
}

# thesis rendered from the stored summaries vs the stored PDF
read_if <- function(f) if (file.exists(rp("verification", f)))
  read.csv(rp("verification", f), check.names = FALSE, na.strings = character()) else NULL
compare <- read_if("render_compare.csv")
numbers <- read_if("render_numbers.csv")
text    <- read_if("render_text.csv")
text_verdicts <- c("identical", "reordered", "URL split elsewhere", "page break elsewhere", "differs")
render_check <- function(display) {
  if (is.null(compare)) return("not compared")
  r <- compare[compare$display == display, ]
  if (!nrow(r)) return("NOT FOUND in the render")
  if (r$identical %in% c(TRUE, "True")) {
    if (r$numbers == 0) "identical (no numbers)" else sprintf("identical (%d numbers)", r$numbers)
  } else sprintf("DIFFERS: %s", gsub("\\|", "vs", r$first_difference))
}

displays$summary_check <- vapply(displays$summaries, summary_check, character(1))
displays$rerun_check   <- vapply(displays$results, rerun_check, character(1))
displays$render_check  <- vapply(displays$display, render_check, character(1))

# the document
md <- c(
  "# Verification",
  "",
  "Written by `replication/verification/make_verification.R` from the logs of the last run of",
  sprintf("`replication/run_checks.R` (%s) on %s, %s, balnet %s and glmnet %s.", ran_at,
          sub("^Running under: ", "", grep("^Running under", session, value = TRUE)),
          sub(" \\(.*", "", session[1]), sub("balnet: ", "", grep("^balnet:", session, value = TRUE)),
          sub("glmnet: ", "", grep("^glmnet:", session, value = TRUE))),
  "",
  "## The checks",
  "",
  "| Link | Check | Result | Source |",
  "|---|---|---|---|",
  sprintf("| Code = the code that ran | `tests/run_tests.R`: the replication code against the code under `R/`, `runs/` and `handoff/`, and against the stored numbers | %d files, %d expectations, %d failures | `logs/tests.txt` |",
          tally[["files"]], tally[["pass"]], tally[["fail"]]),
  if (!is.null(bare))
    sprintf("| Code = the code that ran, without `results/` | `tests/run_tests.R` in a copy of the package without `results/`, the tests reading the rows they use from `tests/fixtures/` | %d files, %d expectations, %d failures, %d skipped (the recomputations in `test_measures.R`, which need `results/`) | `logs/tests_without_results.txt` |",
            bare[["files"]], bare[["pass"]], bare[["fail"]], bare[["skip"]]),
  sprintf("| Raw results = a rerun of the code | `smoke_test.R`: 15 replications of the first and last cell of every grid, n up to 50 000, rerun through `runs/` and compared row by row. Not part of `run_checks.R` (`SMOKE=1` adds it) | %d stored files: %s | `smoke_report.csv` (%s) |",
          nrow(smoke), count(smoke$verdict, verdicts), format(file.mtime(rp("smoke_report.csv")), "%d %b %Y")),
  sprintf("| Summaries = raw results | `verify_summaries.R`: the Study 1 section summaries and the coverage table recomputed from `results/` | %d recomputed, largest difference %s | `logs/summaries.txt` |",
          length(recomputed), format(signif(max(recomputed, na.rm = TRUE), 2))),
  sprintf("| Summaries = raw results | `handoff/scripts/build_summaries.R`: the appendix entries, the Study 2 summaries and `wc_compare` by a second computation | %d of %d passed | `output/thesis_summaries/verification.md` |",
          sum(second == "passed"), length(second)),
  sprintf("| PDF = a render from the summaries | `verify_render.R`: the thesis rendered from the stored summaries in a clean copy. `verification/compare_render.R`: the numbers of every table and figure, every number in the text and the words of every page, against `output/rendered/thesis/AAAAA-Working-Thesis.pdf` | %s | `verification/render_compare.csv`, `verification/render_numbers.csv`, `verification/render_text.csv` |",
          if (is.null(compare) || is.null(text)) "not compared" else
            sprintf("%d of %d displays with identical numbers (%s numbers). In the whole text %d numbers in one PDF and not the other. The %d pages by their words: %s. Both listed below",
                    sum(compare$identical %in% c(TRUE, "True")), nrow(displays),
                    format(sum(compare$numbers), big.mark = " "), sum(numbers$times),
                    nrow(text), count(text$verdict, text_verdicts))),
  "",
  "The steps of the last run:",
  "",
  "```", steps, "```",
  if (any(grepl("^render +failed", steps)) && file.exists(rp("verification", "render_note.md")))
    c("", "The render step stops where there is no TeX. How the render was then made is below."),
  "",
  "## Each table and figure",
  "",
  "Page is the page of the stored PDF. Summary = results gives the largest difference between the",
  "recomputed and the stored numbers, relative to the stored number where that exceeds 1 (above",
  "1e-10 stops `verify_summaries.R`). Rerun = stored counts the `smoke_test.R` verdicts of the results",
  "files behind the display. Render = PDF compares the numbers of the display's rows, in order.",
  "",
  "| Display | Page | qmd chunk | Summary = results | Rerun = stored | Render = PDF |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %s | `%s` (line %s) | %s | %s | %s |", displays$display, displays$page, displays$chunk,
          displays$qmd_line, displays$summary_check, displays$rerun_check, displays$render_check),
  "")

if (file.exists(rp("verification", "render_note.md")))
  md <- c(md, "## How the render was made", "", lines(rp("verification", "render_note.md")), "")
if (!is.null(numbers) && nrow(numbers)) {
  md <- c(md, "## Numbers in one PDF and not the other", "",
          "Counted over the whole text with page numbers dropped. Pages are where the number occurs in that PDF.", "",
          "| Number | Found in | Times | Pages |", "|---|---|---|---|",
          sprintf("| `%s` | %s | %s | %s |", numbers$number, numbers$`in`, numbers$times,
                  ifelse(nchar(numbers$pages) > 40, paste0(substr(numbers$pages, 1, 40), " ..."), numbers$pages)), "")
}
if (!is.null(text)) {
  cell <- function(x) ifelse(nchar(x) == 0, "", paste0("`", ifelse(nchar(x) > 60, paste0(substr(x, 1, 60), " ..."), x), "`"))
  odd <- text[text$verdict != "identical", ]
  md <- c(md, "## Pages whose words are not identical", "",
          "Each page's words read in order, with words split at the end of a line joined. Reordered: the",
          "same words read in another order, as the labels of a figure drawn in another order. URL split",
          "elsewhere: the same letters, a long URL or DOI broken at another point. Page break elsewhere:",
          "the same letters once the next or previous page is added, a line on the other side of a page break.", "",
          "| Page | Verdict | Only in the stored PDF | Only in the render |", "|---|---|---|---|",
          sprintf("| %s | %s | %s | %s |", odd$page, odd$verdict, cell(odd$only_stored), cell(odd$only_render)), "")
}

md <- c(md, "## Not rerun here", "",
        "- The full simulations. The smoke test reruns the first and last cell of every grid.",
        "- The n = 1 million cell (`aipw_misspec_n_dgp2n1mill`), too large for the smoke test.",
        "  `tests/test_runs.R` regenerates its first draw, and the summaries that read it are",
        "  recomputed from the stored file.",
        "")
writeLines(md, rp("VERIFICATION.md"))
cat("written replication/VERIFICATION.md:", nrow(displays), "displays\n")
