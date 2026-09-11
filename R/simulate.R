## R/simulate.R --------------------------------------------------------------

## --- seed streams ----------------------------------------------------------

## L'Ecuyer-CMRG: each stream jumps 2^127 draws ahead, so streams can never
## overlap. Use streams for parrelel reps following Morris, as Default RNG 
## can't guarantee this, Morris (2019)

#' Independent L'Ecuyer-CMRG seed streams
#'
#' Switches the session RNG to L'Ecuyer-CMRG, seeds it, and advances
#' parallel::nextRNGStream() n_streams - 1 times. The stream for a replication
#' is indexed (cell - 1) * num_sim + rep.
#'
#' @param base_seed Integer seed.
#' @param n_streams Number of streams.
#' @return A list of n_streams .Random.seed vectors.
make_seeds <- function(base_seed, n_streams) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(base_seed)
  s <- vector("list", n_streams)
  s[[1]] <- .Random.seed
  for (i in seq_len(n_streams - 1)) s[[i + 1]] <- parallel::nextRNGStream(s[[i]])
  s   # index: (cell - 1) * num_sim + rep
}

## --- one replication --------------------------------------------------------

## Generate, estimate, return the estimate matrix or the error message.

#' One replication: generate a draw and estimate it, catching errors
#'
#' @param cell One-row list of design values passed to dgp_gen.
#' @param dgp_gen Function of a cell returning a draw.
#' @param ... Passed to estimate_all().
#' @return The estimate matrix from estimate_all(), or the error message
#'   (character) if the replication failed.
run_rep <- function(cell, dgp_gen, ...)
  tryCatch(estimate_all(dgp_gen(cell), ...), error = conditionMessage)

## --- driver -----------------------------------------------------------------

#' Run every cell of a design grid in parallel with per-cell checkpoints
#'
#' Sets single-threaded BLAS/OpenMP, opens a future multisession plan and, for
#' each grid row, maps run_rep() over num_sim seed streams with furrr. Each
#' finished cell is written to <stem>_cells/cell_NNN.csv.gz with its stored CV
#' curves in cell_NNN_cv.rds; cells already checkpointed are reloaded and
#' skipped, so an interrupted run resumes. Failed replications keep the
#' estimator block with NA values and the error message in err; a cell whose
#' replications all fail stops the run. At the end the cells are bound and
#' written to out_file, the curves to <stem>_cvcurves.rds, sessionInfo() to
#' <stem>_session.txt, and the checkpoint directory is deleted.
#'
#' @param dgp_gen Function of a cell returning a draw.
#' @param grid Data frame of design cells; must not contain a column "outcome".
#' @param num_sim Replications per cell.
#' @param base_seed Base seed for make_seeds().
#' @param out_file Path of the output csv.gz; the sidecar names derive from it.
#' @param cv_curve_reps Number of leading replications per cell whose full CV
#'   curve is stored.
#' @param workers Number of parallel workers.
#' @param ... Passed to run_rep() and on to estimate_all().
#' @return The long results data.table, invisibly.
simulate_grid <- function(dgp_gen, grid, num_sim, base_seed, out_file,
                          cv_curve_reps = 50,
                          workers = max(1L, future::availableCores(logical = FALSE) - 1L), ...) {
  stopifnot(!"outcome" %in% names(grid))
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1") #prevent oversubscribing threads
  future::plan(future::multisession, workers = workers)  # parallelism at rep level only
  on.exit(future::plan(future::sequential), add = TRUE)
  
  seeds     <- make_seeds(base_seed, nrow(grid) * num_sim)
  res_cells <- vector("list", nrow(grid))
  curves    <- vector("list", nrow(grid))
  
  ck_dir <- sub("\\.csv\\.gz$", "_cells", out_file)   # per-cell checkpoints
  dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)

  for (ci in seq_len(nrow(grid))) {
    ck_file <- file.path(ck_dir, sprintf("cell_%03d.csv.gz", ci))
    cv_file <- file.path(ck_dir, sprintf("cell_%03d_cv.rds", ci))
    if (file.exists(ck_file)) {                     # resume: skip completed cells
      res_cells[[ci]] <- data.table::fread(ck_file, colClasses = list(character = "err"))
      if (file.exists(cv_file)) curves[[ci]] <- readRDS(cv_file)
      next
    }
    t0   <- Sys.time()
    cell <- as.list(grid[ci, , drop = FALSE])
    
    raw <- furrr::future_map(
      seq_len(num_sim),
      function(r) run_rep(cell, dgp_gen, cv_curve = r <= cv_curve_reps, ...),
      .options = furrr::furrr_options(
        seed     = seeds[(ci - 1L) * num_sim + seq_len(num_sim)],
        packages = c("balnet", "glmnet")))
    
    err <- vapply(raw, function(x) if (is.character(x)) x else "", character(1))
    ok  <- err == ""
    if (!any(ok)) stop("cell ", ci, ": all reps failed: ", err[1])
    curves[[ci]] <- lapply(head(raw, cv_curve_reps), attr, "cv_curve")
    
    ## failed reps keep the full estimator block, NA values, err = message
    tmpl <- raw[[which.max(ok)]]
    stopifnot(all(vapply(raw[ok], nrow, integer(1)) == nrow(tmpl)))
    raw[!ok] <- list(tmpl * NA)
    
    res_cells[[ci]] <- data.table::data.table(
      grid[ci, , drop = FALSE],
      sim       = rep(seq_len(num_sim), each = nrow(tmpl)),
      estimator = rep(rownames(tmpl), num_sim),
      do.call(rbind, raw),
      err       = rep(err, each = nrow(tmpl)))
    data.table::fwrite(res_cells[[ci]], ck_file)
    saveRDS(curves[[ci]], cv_file)
    
    cat(sprintf("cell %d/%d: %d failed, %s\n",
                ci, nrow(grid), sum(!ok), format(Sys.time() - t0, digits = 2)))
  }
 
  out <- data.table::rbindlist(res_cells)
  data.table::fwrite(out, out_file)
  saveRDS(curves, sub("\\.csv\\.gz$", "_cvcurves.rds", out_file))
  writeLines(c(capture.output(sessionInfo()),
               paste("balnet:", packageVersion("balnet")),
               paste("glmnet:", packageVersion("glmnet")),
               paste("base_seed:", base_seed),
               paste("date:", Sys.time())),
             sub("\\.csv\\.gz$", "_session.txt", out_file))
  unlink(ck_dir, recursive = TRUE)
  invisible(out)
} 
  
  
  
  
