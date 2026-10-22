# the thesis sources this. attaches sections.R thesis.R thesis_appendix.R as env thesis_api
# its fns look there first so sourced helper files can't overwrite them

local({
  # packages first so they sit below the api on the search path. later ones can't mask them
  for (p in c("dplyr", "tidyr", "purrr", "ggplot2", "gt", "tibble")) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  if ("thesis_api" %in% search()) detach("thesis_api", character.only = TRUE)
  api <- attach(NULL, name = "thesis_api")
  for (f in c("sections.R", "thesis.R", "thesis_appendix.R")) {
    sys.source(here::here("handoff", "R", f), envir = api, keep.source = FALSE)
  }
})
