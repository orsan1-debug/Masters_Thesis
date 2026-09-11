## R/registry.R --------------------------------------------------------------
## One row per completed batch. Design detail lives in the run script; this is
## just the ledger. register_run() is called by the run script after its
## result file exists, never from an on.exit() handler, so a failed or
## interrupted run leaves no row.

registry_cols <- c("batch_id", "component", "dgp", "estimators", "n_sim",
                   "seed", "script", "result_file", "date", "status",
                   "balnet_version", "glmnet_version", "balnet_source")

#' Empty registry ledger
#'
#' @return A zero-row data.frame whose columns are registry_cols, all character.
empty_registry <- function()
  setNames(data.frame(matrix("", 0, length(registry_cols))), registry_cols)

#' Path relative to the repo root (unchanged if outside it)
#'
#' @param p A file path.
#' @return p with the here::here() root prefix removed, forward slashes.
rel_path <- function(p) {
  p    <- normalizePath(p, winslash = "/", mustWork = FALSE)
  root <- paste0(normalizePath(here::here(), winslash = "/"), "/")
  if (startsWith(p, root)) substring(p, nchar(root) + 1) else p
}

#' Install source of the balnet build, if the package metadata exposes one
#'
#' Joins whichever of RemoteSha, RemoteRef, RemoteType and RemoteUrl are
#' present in packageDescription("balnet"); "" when none is.
#'
#' @return A single string.
balnet_source <- function() {
  d <- utils::packageDescription("balnet")
  keep <- c("RemoteSha", "RemoteRef", "RemoteType", "RemoteUrl")
  have <- keep[vapply(keep, function(k) !is.null(d[[k]]), logical(1))]
  if (!length(have)) return("")
  paste(paste0(have, "=", vapply(have, function(k) d[[k]], "")), collapse = "; ")
}

#' Append one completed batch to the registry ledger
#'
#' batch_id is the next row number, date the current time, and the package
#' versions and balnet install source are read from the installed packages;
#' every other field is supplied by the caller. Call it only after the result
#' file has been written.
#'
#' @param component "ipw" or "cv".
#' @param dgp Generator name, for example "dgp2", "gen_data", "dgp_wc_overlap".
#' @param estimators Character vector of estimator names; stored ";"-separated.
#' @param n_sim Replications per cell.
#' @param seed Base or master seed (character if it varies per cell, e.g. "1000 + i").
#' @param script Run script path relative to the repo root.
#' @param result_file Result file(s) relative to the repo root; a glob for
#'   per-cell batches.
#' @param status Free text; "complete" by default.
#' @param path Path of the registry csv.
#' @return The batch_id, invisibly.
register_run <- function(component, dgp, estimators, n_sim, seed, script,
                         result_file, status = "complete",
                         path = here::here("registry.csv")) {
  reg <- if (file.exists(path)) read.csv(path, colClasses = "character")
         else empty_registry()
  entry <- data.frame(batch_id       = as.character(nrow(reg) + 1L),
                      component      = component,
                      dgp            = dgp,
                      estimators     = paste(estimators, collapse = ";"),
                      n_sim          = as.character(n_sim),
                      seed           = as.character(seed),
                      script         = script,
                      result_file    = rel_path(result_file),
                      date           = format(Sys.time(), "%Y-%m-%d %H:%M"),
                      status         = status,
                      balnet_version = as.character(utils::packageVersion("balnet")),
                      glmnet_version = as.character(utils::packageVersion("glmnet")),
                      balnet_source  = balnet_source(),
                      stringsAsFactors = FALSE)
  write.csv(rbind(reg, entry), path, row.names = FALSE)
  invisible(entry$batch_id)
}

#' Run a simulation grid and, on success, log it in the registry
#'
#' Wraps simulate_grid(). register_run() runs only after simulate_grid() has
#' returned, so a failed run leaves no registry row. The estimator list is
#' read off the result (diagnostic rows excluded); status is "complete",
#' followed by the count of replications whose err column is non-empty when
#' there are any, and by note when given.
#'
#' @param dgp_gen Function of one grid cell returning a draw.
#' @param grid Data frame of design cells.
#' @param num_sim Replications per cell.
#' @param base_seed Base seed for the L'Ecuyer-CMRG streams.
#' @param out_file Output csv.gz path passed to simulate_grid().
#' @param dgp Generator name recorded in the registry.
#' @param script Run script path recorded in the registry.
#' @param note Optional text appended to status (for example "patched cv.balnet").
#' @param registry Path of the registry csv.
#' @return The long results data.table from simulate_grid().
run_batch <- function(dgp_gen, grid, num_sim, base_seed, out_file,
                      dgp, script, note = NULL,
                      registry = here::here("registry.csv")) {
  res   <- simulate_grid(dgp_gen, grid, num_sim, base_seed, out_file)
  n_err <- sum(res$err != "")
  ests  <- unique(res$estimator)
  ests  <- ests[!grepl("^(lam_|nnz_|smd[01]_|cvloss_|trunc05|prev|emin|emax|nout)", ests)]
  status <- paste(c("complete",
                    if (n_err > 0) paste0("err rows: ", n_err),
                    note), collapse = "; ")
  register_run(component = "ipw", dgp = dgp, estimators = ests,
               n_sim = num_sim, seed = base_seed, script = script,
               result_file = out_file, status = status, path = registry)
  res
}
