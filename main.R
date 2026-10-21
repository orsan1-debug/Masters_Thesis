# master file. Rscript main.R from the project root (or source it in RStudio)
# set the options below. defaults need no results/ and take about 10 min

reps      <- 15      # reps per cell for reruns. NA = the reps we ran (printed below)
max_n     <- 50000   # rerun cells up to this n. Inf = all incl n = 1 million
all_cells <- FALSE   # FALSE = first and last cell of each grid only
workers   <- max(1L, parallel::detectCores() - 1L)
rerun     <- FALSE   # TRUE reruns the sims. at NA reps this is weeks of machine time
render    <- TRUE    # needs Quarto XeLaTeX and Times New Roman

rscript <- function(script, env = character()) {
  status <- withr::with_envvar(env, system2(file.path(R.home("bin"), "Rscript"),
                                            shQuote(here::here(script))))
  if (status != 0) stop(script, " stopped with status ", status)
}
has_results <- dir.exists(here::here("results", "aipw"))

rscript("replication/00_setup.R")

# reps we ran per grid. read from the run scripts
grids_of <- function(f) local({ source(here::here("replication", "runs", f), local = TRUE); grids })
s1 <- sapply(grids_of("study1.R"), `[[`, "replications")
s2 <- sapply(grids_of("study2.R"), \(g) if (is.null(g$replications)) 2500 else g$replications)
print(data.frame(study = c(rep("study 1", length(s1)), rep("study 2", length(s2)), "wong chan"),
                 grid = c(names(s1), names(s2), "wong_chan"), reps = c(s1, s2, 1000)),
      row.names = FALSE)

rscript("replication/tests/run_tests.R")
if (has_results) rscript("replication/verify_summaries.R")
if (render) rscript("replication/verify_render.R")   # stored pdf left alone

# reruns go to output/rerun/ then get compared row by row w/ results/ if its there
if (rerun) {
  out   <- here::here("output", "rerun", if (is.na(reps)) "reps_full" else paste0("reps_", reps))
  env   <- c(WORKERS = workers, MAX_N = max_n, CELLS = if (all_cells) "all" else "ends",
             N_SIM = if (is.na(reps)) NULL else reps)
  parts <- c(study1 = "study1.R", study2 = "study2.R", wc = "wong_chan.R")
  for (p in names(parts))
    rscript(file.path("replication", "runs", parts[[p]]), c(env, OUT_DIR = file.path(out, p)))
  if (has_results)
    rscript("replication/smoke_test.R",
            c(SKIP_RUNS = "1", OUT_DIR = out, REPORT = file.path(out, "compare.csv")))
}
