# runs a design grid in parallel w/ per cell checkpoints
# one L'Ecuyer-CMRG stream per rep. 2^127 apart so reps never overlap (Morris et al 2019 s5)

#' independent L'Ecuyer-CMRG seed streams. list of .Random.seed vectors
make_seeds <- function(base_seed, n_streams) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(base_seed)
  s <- vector("list", n_streams)
  s[[1]] <- .Random.seed
  for (i in seq_len(n_streams - 1)) s[[i + 1]] <- parallel::nextRNGStream(s[[i]])
  s
}

#' one rep. draw + estimate. returns the error message if it fails
run_rep <- function(cell, dgp_gen, estimator = estimate_all, ...)
  tryCatch(estimator(dgp_gen(cell), ...), error = conditionMessage)

#' run every grid cell in parallel. each cell saved to <stem>_cells/ so an interrupted run resumes
#' failed reps get NA + message in err. a cell w/ all reps failed stops the run
simulate_grid <- function(dgp_gen, grid, num_sim, base_seed, out_file,
                          estimator = estimate_all,  # add-only run = other fn + new out_file
                          workers = max(1L, future::availableCores(logical = FALSE) - 1L), ...) {
  stopifnot(!"outcome" %in% names(grid))
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")  # no oversubscription
  strategy <- getOption("sim.plan", future::multisession)  # eg future::multicore on macOS
  future::plan(strategy, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  seeds     <- make_seeds(base_seed, nrow(grid) * num_sim)
  res_cells <- vector("list", nrow(grid))
  ck_dir    <- sub("\\.csv\\.gz$", "_cells", out_file)
  dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (ci in seq_len(nrow(grid))) {
    ck_file <- file.path(ck_dir, sprintf("cell_%03d.csv.gz", ci))
    if (file.exists(ck_file)) {                     # resume. skip done cells
      res_cells[[ci]] <- data.table::fread(ck_file, colClasses = list(character = "err"))
      next
    }
    t0   <- Sys.time()
    cell <- as.list(grid[ci, , drop = FALSE])
    
    raw <- furrr::future_map(
      seq_len(num_sim),
      function(r) run_rep(cell, dgp_gen, estimator, ...),
      .options = furrr::furrr_options(
        seed     = seeds[(ci - 1L) * num_sim + seq_len(num_sim)],
        packages = c("balnet", "glmnet")))
    
    err <- vapply(raw, function(x) if (is.character(x)) x else "", character(1))
    ok  <- err == ""
    if (!any(ok)) stop("cell ", ci, ": all reps failed: ", err[1])
    
    tmpl <- raw[[which.max(ok)]]                    # failed reps get an NA block + err
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
  
  
  
  
