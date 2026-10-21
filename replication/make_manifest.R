# lists the thesis files their deps and the test fixtures w/ size and md5 so a copy can be checked
#   Rscript replication/make_manifest.R   writes replication/manifest.csv

roles <- list(
  thesis       = c("thesis/AAAAA Working Thesis.qmd", "thesis/references.bib", "_quarto.yml",
                   "output/rendered/thesis/AAAAA-Working-Thesis.pdf"),
  replication  = list.files(here::here("replication"), recursive = TRUE, full.names = FALSE,
                            pattern = "\\.(R|md)$") |> file.path("replication", . = _),
  record       = c(list.files(here::here("R"), "\\.R$", full.names = FALSE) |> file.path("R", . = _),
                   list.files(here::here("runs"), "\\.R$", recursive = TRUE, full.names = FALSE) |>
                     file.path("runs", . = _),
                   list.files(here::here("handoff"), "\\.R$", recursive = TRUE, full.names = FALSE) |>
                     file.path("handoff", . = _)),
  package      = "pkg/balnet_0.0.4.tar.gz",
  fixtures     = list.files(here::here("replication", "tests", "fixtures"), recursive = TRUE,
                            full.names = FALSE) |> file.path("replication", "tests", "fixtures", . = _),
  summaries    = list.files(here::here("output", "thesis_summaries"), full.names = FALSE) |>
                   file.path("output", "thesis_summaries", . = _),
  results      = list.files(here::here("results"), "\\.(csv\\.gz|rds|txt|csv)$", recursive = TRUE,
                            full.names = FALSE) |> file.path("results", . = _))

# record code for designs the thesis does not use stays with the author
not_thesis <- c("^runs/ipw/", "^runs/cv/basic_dgp/", "^runs/cv/wz_replication/run_wc_(basis|overlap|overlap_n|path)\\.R$",
                "^handoff/style_examples\\.R$", "^handoff/scripts/(compare_legacy|tidy_copy)\\.R$")
roles$record <- grep(paste(not_thesis, collapse = "|"), roles$record, value = TRUE, invert = TRUE)

manifest <- do.call(rbind, lapply(names(roles), function(role) {
  files <- roles[[role]]
  files <- files[file.exists(here::here(files))]
  if (!length(files)) return(NULL)
  data.frame(role = role, path = files,
             bytes = file.size(here::here(files)),
             md5 = unname(tools::md5sum(here::here(files))),
             modified = format(file.mtime(here::here(files)), "%Y-%m-%d %H:%M"))
}))
write.csv(manifest, here::here("replication", "manifest.csv"), row.names = FALSE)
cat(nrow(manifest), "files,", round(sum(manifest$bytes) / 1e6), "MB. Written replication/manifest.csv\n")
print(aggregate(bytes ~ role, manifest, function(b) c(files = length(b), MB = round(sum(b) / 1e6))))
