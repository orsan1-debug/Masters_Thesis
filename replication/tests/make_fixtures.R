# cuts the rows the tests use from results/ into tests/fixtures/ so tests run w/o the 7 GB results/
#   Rscript replication/tests/make_fixtures.R   from the project root. needs results/

library(data.table)

out <- here::here("replication", "tests", "fixtures")
study1 <- local({ source(here::here("replication", "runs", "study1.R"), local = TRUE); grids })
study2 <- local({ source(here::here("replication", "runs", "study2.R"), local = TRUE); grids })
stored <- function(g) unname(unlist(g$stored))

cut_rows <- function(from, keep, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  d <- fread(from, colClasses = list(character = "err"))
  fwrite(d[keep(d)], to)
}

# study 1. oracle_hajek rep 1 of every stored file + a few reps of the sample size files
sample_size <- c("aipw_correct_n_dgp2", "aipw_add_correct_n400k_dgp2",
                 "aipw_misspec_n_dgp2", "aipw_add_misspec_n_dgp2")
for (stem in unique(unlist(lapply(study1, stored)))) {
  keep <- function(d) {
    rows <- d$estimator == "oracle_hajek" & d$sim == 1
    if ("batch" %in% names(d)) rows <- rows & d$batch == 1
    if (stem %in% sample_size) rows <- rows | (d$sim <= 2 & d$n <= 10000)
    if (stem == "aipw_misspec_n_dgp2") rows <- rows | (d$n == 500 & d$sim <= 3)
    rows
  }
  cut_rows(here::here("results", "aipw", paste0(stem, ".csv.gz")), keep,
           file.path(out, "aipw", paste0(stem, ".csv.gz")))
}

# study 2. lambda grids of overlap folds and denser + rep 1 of the folds and denser files
r2500 <- function(root, stem, ext) file.path(root, "cv", "our dgp", "r2500", paste0(stem, ext))
for (name in c("overlap", "folds", "denser")) for (stem in stored(study2[[name]])) {
  dir.create(dirname(r2500(out, stem, "")), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(r2500(here::here("results"), stem, "_lamgrid.csv"), r2500(out, stem, "_lamgrid.csv"),
                      overwrite = TRUE))
  if (name != "overlap")
    cut_rows(r2500(here::here("results"), stem, ".csv.gz"), function(d) d$sim == 1, r2500(out, stem, ".csv.gz"))
}

# wong chan. the cell sigma = 10 c = 2 w/ its first rep
wc <- readRDS(here::here("results", "cv", "wz_replication", "wc_noise_s10_c2.rds"))
wc$res <- wc$res[1]
dir.create(file.path(out, "cv", "wz_replication"), recursive = TRUE, showWarnings = FALSE)
saveRDS(wc, file.path(out, "cv", "wz_replication", "wc_noise_s10_c2.rds"))

files <- list.files(out, recursive = TRUE)
cat(length(files), "files,", round(sum(file.size(file.path(out, files))) / 1e6, 1),
    "MB written to replication/tests/fixtures\n")
