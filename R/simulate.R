## R/simulate.R --------------------------------------------------------------

## --- seed streams ----------------------------------------------------------

## L'Ecuyer-CMRG streams are 2^127 draws apart, so parallel reps never overlap
## (Morris et al. 2019 s.5, seeds). Stream index: (cell - 1) * num_sim + rep.

#' Independent L'Ecuyer-CMRG seed streams
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
  s
}

## --- one replication --------------------------------------------------------

#' One replication: generate a draw, estimate it, catch errors
#'
#' @param cell One-row list of design values passed to dgp_gen.
#' @param dgp_gen Function of a cell returning a draw.
#' @param ... Passed to estimate_all().
#' @return The estimate matrix, or the error message (character) on failure.
run_rep <- function(cell, dgp_gen, ...)
  tryCatch(estimate_all(dgp_gen(cell), ...), error = conditionMessage)

## --- driver -----------------------------------------------------------------

#' Run every cell of a design grid in parallel with per-cell checkpoints
#'
#' Single-threaded BLAS/OpenMP, one future multisession plan; each grid row is
#' mapped over num_sim seed streams with furrr and written to
#' <stem>_cells/cell_NNN.csv.gz, so an interrupted run resumes. Failed reps
#' keep the estimator block with NA values and the message in err; a cell whose
#' reps all fail stops the run. At the end the cells are bound into out_file,
#' sessionInfo() goes to <stem>_session.txt and the checkpoints are deleted.
#'
#' @param dgp_gen Function of a cell returning a draw.
#' @param grid Data frame of design cells; must not contain a column "outcome".
#' @param num_sim Replications per cell.
#' @param base_seed Base seed for make_seeds().
#' @param out_file Path of the output csv.gz; sidecar names derive from it.
#' @param workers Number of parallel workers.
#' @param ... Passed to run_rep() and on to estimate_all().
#' @return The long results data.table, invisibly.
simulate_grid <- function(dgp_gen, grid, num_sim, base_seed, out_file,
                          workers = max(1L, future::availableCores(logical = FALSE) - 1L), ...) {
  stopifnot(!"outcome" %in% names(grid))
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")  # no oversubscription
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  seeds     <- make_seeds(base_seed, nrow(grid) * num_sim)
  res_cells <- vector("list", nrow(grid))
  ck_dir    <- sub("\\.csv\\.gz$", "_cells", out_file)
  dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (ci in seq_len(nrow(grid))) {
    ck_file <- file.path(ck_dir, sprintf("cell_%03d.csv.gz", ci))
    if (file.exists(ck_file)) {                     # resume: skip completed cells
      res_cells[[ci]] <- data.table::fread(ck_file, colClasses = list(character = "err"))
      next
    }
    t0   <- Sys.time()
    cell <- as.list(grid[ci, , drop = FALSE])
    
    raw <- furrr::future_map(
      seq_len(num_sim),
      function(r) run_rep(cell, dgp_gen, ...),
      .options = furrr::furrr_options(
        seed     = seeds[(ci - 1L) * num_sim + seq_len(num_sim)],
        packages = c("balnet", "glmnet")))
    
    err <- vapply(raw, function(x) if (is.character(x)) x else "", character(1))
    ok  <- err == ""
    if (!any(ok)) stop("cell ", ci, ": all reps failed: ", err[1])
    
    tmpl <- raw[[which.max(ok)]]                    # failed reps: NA block + err
    stopifnot(all(vapply(raw[ok], nrow, integer(1)) == nrow(tmpl)))
    raw[!ok] <- list(tmpl * NA)
    
    res_cells[[ci]] <- data.table::data.table(
      grid[ci, , drop = FALSE],
      sim       = rep(seq_len(num_sim), each = nrow(tmpl)),
      estimator = rep(rownames(tmpl), num_sim),
      do.call(rbind, raw),
      err       = rep(err, each = nrow(tmpl)))
    data.table::fwrite(res_cells[[ci]], ck_file)
    
    cat(sprintf("cell %d/%d: %d failed, %s\n",
                ci, nrow(grid), sum(!ok), format(Sys.time() - t0, digits = 2)))
  }
  
  out <- data.table::rbindlist(res_cells)
  data.table::fwrite(out, out_file)
  writeLines(c(capture.output(sessionInfo()),
               paste("balnet:", packageVersion("balnet")),
               paste("glmnet:", packageVersion("glmnet")),
               paste("base_seed:", base_seed),
               paste("date:", Sys.time())),
             sub("\\.csv\\.gz$", "_session.txt", out_file))
  unlink(ck_dir, recursive = TRUE)
  invisible(out)
}
  
  
  
  
