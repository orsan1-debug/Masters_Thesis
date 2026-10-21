# ledger of completed batches. one row each. design detail lives in the run script
# call register_run() after the result file exists not in on.exit() so failed runs leave no row

registry_cols <- c("batch_id", "component", "dgp", "estimators", "n_sim",
                   "seed", "script", "result_file", "date", "status",
                   "balnet_version", "glmnet_version", "balnet_source")

#' empty ledger. zero rows and all cols character
empty_registry <- function()
  setNames(data.frame(matrix("", 0, length(registry_cols))), registry_cols)

#' path relative to the repo root (unchanged if outside it)
rel_path <- function(p) {
  p    <- normalizePath(p, winslash = "/", mustWork = FALSE)
  root <- paste0(normalizePath(here::here(), winslash = "/"), "/")
  if (startsWith(p, root)) substring(p, nchar(root) + 1) else p
}

#' balnet install source from the Remote* fields of its packageDescription. "" if none
balnet_source <- function() {
  d <- utils::packageDescription("balnet")
  keep <- c("RemoteSha", "RemoteRef", "RemoteType", "RemoteUrl")
  have <- keep[vapply(keep, function(k) !is.null(d[[k]]), logical(1))]
  if (!length(have)) return("")
  paste(paste0(have, "=", vapply(have, function(k) d[[k]], "")), collapse = "; ")
}

#' append one completed batch. batch_id = next row number and pkg versions read from the install
#' seed can be text if it varies per cell (eg "1000 + i"). result_file can be a glob
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
