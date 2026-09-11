## R/registry.R --------------------------------------------------------------
## One row per run, appended automatically by run_batch().
## Design detail lives in the run script; this is just the ledger.

registry_cols <- c("batch_id", "batch_file", "label", "seed",
                   "balnet_version", "date_run", "status",
                   "superseded_by", "verified", "notes")

empty_registry <- function()
  setNames(data.frame(matrix("", 0, length(registry_cols))), registry_cols)

register_run <- function(entry, path = here::here("registry.csv")) {
  reg <- if (file.exists(path)) read.csv(path, colClasses = "character")
  else empty_registry()
  
  entry <- modifyList(list(
    batch_id       = nrow(reg) + 1L,
    date_run       = format(Sys.time(), "%Y-%m-%d %H:%M"),
    balnet_version = as.character(packageVersion("balnet")),
    verified       = "no"), entry)
  
  entry[setdiff(registry_cols, names(entry))] <- ""
  write.csv(rbind(reg, as.data.frame(entry[registry_cols])),
            path, row.names = FALSE)
  invisible(entry$batch_id)
}

run_batch <- function(dgp_gen, grid, num_sim, base_seed, out_file,
                      meta = list(), registry = here::here("registry.csv")) {
  status <- "failed"; n_err <- NA
  
  on.exit(register_run(c(meta, list(
    batch_file = basename(out_file), seed = base_seed,
    status = status, notes = paste0("err rows: ", n_err))), registry),
    add = TRUE)
  
  res    <- simulate_grid(dgp_gen, grid, num_sim, base_seed, out_file)
  n_err  <- sum(res$err != "")   # ADAPT: your err column
  status <- "complete"
  res
}