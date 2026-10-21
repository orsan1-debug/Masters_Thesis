# study2_rules() gives the verdict of s2_cell_flags() (handoff/R/thesis.R) in every stored cell

library(testthat)
source(here::here("replication", "R", "measures.R"))

# only the 2 objects needed from handoff/R/thesis.R so gt is not required
record <- new.env()
library(dplyr); library(tidyr)
for (e in parse(here::here("handoff", "R", "thesis.R"))) {
  target <- if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=")) as.character(e[[2]]) else ""
  if (target %in% c("s2_cell_flags", "s2_selectors")) eval(e, envir = record)
}

summaries <- list.files(here::here("output", "thesis_summaries"), "^cv_.*\\.rds$",
                        full.names = TRUE)

test_that("the two rules agree with the thesis code in every stored Study 2 cell", {
  if (length(summaries) == 0) skip("no Study 2 summary present")
  for (f in summaries) {
    s <- readRDS(f)
    ours   <- study2_rules(s$perf, keys = c("part", "cell", "outcome", "strength", "est"))
    theirs <- record$s2_cell_flags(s)
    j <- inner_join(ours, theirs, by = c("part", "cell", "outcome", "strength", "est"),
                    suffix = c("", ".thesis"))
    expect_equal(nrow(j), nrow(theirs), label = basename(f))
    expect_identical(j$interior, j$interior.thesis, label = paste(basename(f), "interior"))
    expect_identical(j$win, j$win.thesis, label = paste(basename(f), "win"))
    expect_identical(j$lead, j$lead.thesis, label = paste(basename(f), "lead"))
    expect_equal(j$path_min, j$path_min.thesis, tolerance = 1e-12)
    expect_equal(j$end, j$end.thesis, tolerance = 1e-12)
  }
})

test_that("the rules behave as defined on made-up cells", {
  cell <- function(arm, rmse, mcse, k = NA) tibble(spec = "s", level = "a", outcome = "linear",
                                                   strength = 1, est = "ipw", arm = arm,
                                                   k = k, reps = 2500, rmse = rmse, mcse_rmse = mcse)
  perf <- bind_rows(cell("fixed", 1.00, 0.01),
                    cell("path", 0.96, 0.01, 1), cell("path", 1.00, 0.01, 2),
                    cell("cv.bloss", 0.97, 0.01), cell("boot.inf", 0.99, 0.01))
  r <- study2_rules(perf)
  expect_false(r$interior)                 # 4% below is not 5%
  expect_false(r$win)                      # 0.03 below is not 3 sqrt(2) 0.01 = 0.042
  expect_equal(r$lead, "cv.bloss")
  perf2 <- perf |> mutate(rmse = ifelse(arm == "path" & k == 1, 0.94, rmse),
                          rmse = ifelse(arm == "cv.bloss", 0.95, rmse))
  r2 <- study2_rules(perf2)
  expect_true(r2$interior)
  expect_true(r2$win)
})
