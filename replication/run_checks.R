# runs each check in its own R process. logs to replication/logs/. smoke test ~20 min on 17 workers
#   Rscript replication/run_checks.R   about 13 minutes and needs results/. SMOKE=1 adds smoke test

steps <- c(setup = "00_setup.R", tests = "tests/run_tests.R",
           summaries = "verify_summaries.R", render = "verify_render.R",
           manifest = "make_manifest.R")
if (Sys.getenv("SMOKE") == "1") steps <- c(steps, smoke = "smoke_test.R")
logs <- here::here("replication", "logs")
dir.create(logs, showWarnings = FALSE)

status <- vapply(names(steps), function(step) {
  started <- Sys.time()
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                  shQuote(here::here("replication", steps[[step]])),
                                  stdout = TRUE, stderr = TRUE))
  code <- if (is.null(attr(out, "status"))) 0 else attr(out, "status")
  writeLines(out, file.path(logs, paste0(step, ".txt")))
  cat(sprintf("%-10s %-20s %s\n", step, if (code == 0) "ok" else paste("failed, status", code),
              format(Sys.time() - started, digits = 2)))
  code
}, numeric(1))
if (any(status != 0)) stop("see replication/logs/ for the failed steps")
