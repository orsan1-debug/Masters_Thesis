# master file. Rscript main.R from the project root (or source it in RStudio)
# set the options below. defaults need no results/ and take about 10 min

reps      <- 15      # reps per cell for reruns. NA = the reps we ran (printed below)
max_n     <- 50000   # rerun cells up to this n. Inf = all incl n = 1 million
all_cells <- FALSE   # FALSE = first and last cell of each grid only
workers   <- max(1L, parallel::detectCores() - 1L)   # at most 6 for n >= 100k w/ 24 GB
rerun     <- FALSE   # TRUE reruns the sims. at NA reps this is weeks of machine time
render    <- TRUE    # needs Quarto XeLaTeX and Times New Roman

if (!file.exists(".here")) stop("run main.R from the project root")
# base R only as 00_setup.R installs the packages. NA in env unsets that var
rscript <- function(script, env = character()) {
  if (any(!is.na(env))) do.call(Sys.setenv, as.list(env[!is.na(env)]))
  Sys.unsetenv(names(env)[is.na(env)])
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script))
  if (status != 0) stop(script, " stopped with status ", status)
}
out <- file.path(getwd(), "output", "rerun")   # all main.R writes goes here
dir.create(out, recursive = TRUE, showWarnings = FALSE)
has_results <- dir.exists(file.path("results", "aipw"))

rscript("replication/00_setup.R", c(OUT_DIR = out))

# reps we ran per grid. the study 2 default 2500 and wong chan 1000 are as in their run scripts
grids_of <- function(f) local({ source(file.path("replication", "runs", f), local = TRUE); grids })
s1 <- sapply(grids_of("study1.R"), `[[`, "replications")
s2 <- sapply(grids_of("study2.R"), \(g) if (is.null(g$replications)) 2500 else g$replications)
print(data.frame(study = c(rep("study 1", length(s1)), rep("study 2", length(s2)), "wong chan"),
                 grid = c(names(s1), names(s2), "wong_chan"), reps = c(s1, s2, 1000)),
      row.names = FALSE)

rscript("replication/tests/run_tests.R")
if (has_results) rscript("replication/verify_summaries.R")

# reruns then a row by row compare w/ results/ if its there. one folder per setting
if (rerun) {
  runs  <- file.path(out, paste0("reps_", if (is.na(reps)) "full" else reps, "_n", max_n,
                                 if (all_cells) "_all" else "_ends"))
  env   <- c(WORKERS = workers, MAX_N = max_n, CELLS = if (all_cells) "all" else "ends",
             N_SIM = reps, GRIDS = NA)
  parts <- c(study1 = "study1.R", study2 = "study2.R", wc = "wong_chan.R")
  for (p in names(parts))
    rscript(file.path("replication", "runs", parts[[p]]), c(env, OUT_DIR = file.path(runs, p)))
  if (has_results)
    rscript("replication/smoke_test.R",
            c(SKIP_RUNS = "1", PARTS = NA, OUT_DIR = runs, REPORT = file.path(runs, "compare.csv")))
}

# last as it needs TeX. the stored pdf and verification/ are left alone
if (render) rscript("replication/verify_render.R", c(OUT_DIR = file.path(out, "render")))
