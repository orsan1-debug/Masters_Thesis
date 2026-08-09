## R/simulate.R --------------------------------------------------------------

## --- seed streams ----------------------------------------------------------

## L'Ecuyer-CMRG: each stream jumps 2^127 draws ahead, so streams can never
## overlap. Use streams for parrelel reps following Morris, as Default RNG 
## can't guarantee this, Morris (2019)

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

run_rep <- function(cell, dgp_gen, ...)
  tryCatch(estimate_all(dgp_gen(cell), ...), error = conditionMessage)

## --- driver -----------------------------------------------------------------

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
  
  
  
  