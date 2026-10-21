# runs all test files here. w/o results/ they use tests/fixtures/ and 3 test_measures.R checks skip
#   Rscript replication/tests/run_tests.R   from the project root. writes only to tempdir

library(testthat)
tests <- sort(list.files(here::here("replication", "tests"), "^test_.*\\.R$",
                         full.names = TRUE))
results <- lapply(tests, function(f) {
  cat("\n==", basename(f), "==\n")
  test_file(f, reporter = "progress")
})
failed <- sum(vapply(results, function(r) { d <- as.data.frame(r); sum(d$failed) + sum(d$error) }, numeric(1)))
cat(sprintf("\n%d test file(s), %d failure(s) or error(s)\n", length(tests), failed))
if (failed > 0) stop("some checks failed")
