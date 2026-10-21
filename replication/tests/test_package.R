# the package is whole and its indexes agree (manifest.csv run scripts and displays.csv)
# runs with or without results/

library(testthat)
manifest <- read.csv(here::here("replication", "manifest.csv"))
displays <- read.csv(here::here("replication", "verification", "displays.csv"))

test_that("the thesis, the record code, balnet, the summaries and the fixtures are the manifest's", {
  m <- manifest[manifest$role %in% c("thesis", "record", "package", "summaries", "fixtures"), ]
  paths <- here::here(m$path)
  here_ <- file.exists(paths)
  expect_true(all(here_), label = paste("missing:", paste(m$path[!here_], collapse = ", ")))
  expect_equal(unname(file.size(paths[here_])), m$bytes[here_])
  expect_equal(unname(tools::md5sum(paths[here_])), m$md5[here_])
})

test_that("every stored file the run scripts list is in the manifest, with its session", {
  listed <- function(script) {
    grids <- local({ source(here::here("replication", "runs", script), local = TRUE); grids })
    unique(unname(unlist(lapply(grids, function(g) g$stored))))
  }
  in_manifest <- function(paths) paths %in% manifest$path
  study1 <- file.path("results", "aipw", listed("study1.R"))
  study2 <- listed("study2.R")
  r2500  <- file.path("results", "cv", "our dgp", "r2500", study2)
  r1000  <- file.path("results", "cv", "our dgp", "r1000", study2)
  study2 <- ifelse(in_manifest(paste0(r2500, ".csv.gz")), r2500, r1000)
  for (ext in c(".csv.gz", "_session.txt"))
    expect_true(all(in_manifest(paste0(c(study1, study2), ext))), label = paste("stored", ext))
  expect_true(all(in_manifest(paste0(study2, "_lamgrid.csv"))), label = "Study 2 lambda grids")
  wong_chan <- sprintf("results/cv/wz_replication/wc_noise_s%g_c%g.rds",
                       rep(c(1, 10, 30, 100), 4), rep(1:4, each = 4))
  expect_true(all(in_manifest(wong_chan)), label = "Wong and Chan cells")
})

test_that("every display is the labelled cell at its line of the thesis qmd", {
  qmd <- readLines(here::here("thesis", "AAAAA Working Thesis.qmd"), warn = FALSE)
  at  <- qmd[displays$qmd_line]
  ok  <- at == paste("#| label:", displays$chunk) |
         startsWith(at, paste0("::: {#", displays$chunk, "}")) |
         startsWith(at, paste0("::: {#", displays$chunk, " "))
  ok  <- ok %in% TRUE
  expect_true(all(ok), label = paste("not at their line:", paste(displays$display[!ok], collapse = ", ")))
})

test_that("every stored summary a display reads is here and in the manifest", {
  ids <- displays$summaries[!startsWith(displays$summaries, "design")]
  ids <- unique(sub(":.*", "", trimws(unlist(strsplit(ids, ";")))))
  files <- file.path("output", "thesis_summaries", paste0(ids, ".rds"))
  expect_true(all(file.exists(here::here(files))), label = "summaries here")
  expect_true(all(files %in% manifest$path), label = "summaries in the manifest")
})
