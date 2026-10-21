# grid loop on L'Ecuyer-CMRG streams (Methods 3.2). rep r in block b uses stream (b - 1) * R + r
# R = the grid's reps. cells in a block share draws and a short run gets the first reps

#' first n L'Ecuyer-CMRG streams from a base seed
seed_streams <- function(base_seed, n) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(base_seed)
  streams <- vector("list", n)
  streams[[1]] <- .Random.seed
  for (i in seq_len(n - 1)) streams[[i + 1]] <- parallel::nextRNGStream(streams[[i]])
  streams
}

#' one rep from one stream. draw then estimate. returns the error message on failure
run_replication <- function(stream, cell, generate, estimate) {
  assign(".Random.seed", stream, envir = globalenv())
  tryCatch(estimate(generate(cell)), error = conditionMessage)
}

#' every cell of a grid into one results file. failed reps kept as NA rows w/ the error
#' cells are written as they finish so a stopped run resumes. stream_block = the full run's reps
run_grid <- function(generate, grid, replications, base_seed, out_file, estimate,
                     workers = 1, stream_block = replications, cells = seq_len(nrow(grid))) {
  stopifnot(!"outcome" %in% names(grid), replications <= stream_block)
  if (!length(cells)) return(invisible(NULL))
  block    <- stream_blocks(grid)
  design   <- grid[setdiff(names(grid), "block")]
  streams  <- seed_streams(base_seed, max(block[cells]) * stream_block)
  cell_dir <- sub("\\.csv\\.gz$", "_cells", out_file)
  dir.create(cell_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")   # one thread per worker

  cluster <- NULL
  if (workers > 1) {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(root)
      for (f in c("dgp.R", "estimators.R", "selectors.R")) source(file.path(root, "replication", "R", f)),
      root = here::here())
    parallel::clusterExport(cluster, c("generate", "estimate", "run_replication"), envir = environment())
  }
  over_reps <- function(X, ...) if (is.null(cluster)) lapply(X, ...) else parallel::parLapply(cluster, X, ...)

  done <- vector("list", nrow(grid))
  for (i in cells) {
    cell_file <- file.path(cell_dir, sprintf("cell_%03d.csv.gz", i))
    if (file.exists(cell_file)) {
      done[[i]] <- data.table::fread(cell_file, colClasses = list(character = "err"))
      next
    }
    started <- Sys.time()
    at  <- (block[i] - 1) * stream_block + seq_len(replications)
    raw <- over_reps(streams[at], run_replication, cell = as.list(design[i, , drop = FALSE]),
                     generate = generate, estimate = estimate)
    err <- vapply(raw, function(x) if (is.character(x)) x else "", character(1))
    ok  <- err == ""
    if (!any(ok)) stop("cell ", i, ": every replication failed: ", err[1])
    template <- raw[[which.max(ok)]]
    raw[!ok] <- list(template * NA)
    done[[i]] <- data.table::data.table(
      design[i, , drop = FALSE],
      sim       = rep(seq_len(replications), each = nrow(template)),
      estimator = rep(rownames(template), replications),
      do.call(rbind, raw),
      err       = rep(err, each = nrow(template)))
    data.table::fwrite(done[[i]], cell_file)
    cat(sprintf("cell %d of %d: %d failed, %s\n", i, nrow(grid), sum(!ok),
                format(Sys.time() - started, digits = 2)))
  }

  out <- data.table::rbindlist(done)
  data.table::fwrite(out, out_file)
  writeLines(c(capture.output(sessionInfo()),
               paste("balnet:", packageVersion("balnet")), paste("glmnet:", packageVersion("glmnet")),
               paste("base_seed:", base_seed), paste("date:", Sys.time())),
             sub("\\.csv\\.gz$", "_session.txt", out_file))
  unlink(cell_dir, recursive = TRUE)
  invisible(out)
}

#' a grid's stream blocks. the block column or else the row numbers
stream_blocks <- function(grid) if (is.null(grid$block)) seq_len(nrow(grid)) else grid$block

#' each design point under each value given. one cell each
#' cells of a point share its stream block (its row unless there is a block column) so same draws
under_each <- function(points, ...) {
  values <- expand.grid(..., stringsAsFactors = FALSE)
  if (is.null(points$block)) points$block <- seq_len(nrow(points))
  cells <- cbind(points[rep(seq_len(nrow(points)), each = nrow(values)), , drop = FALSE],
                 values[rep(seq_len(nrow(values)), nrow(points)), , drop = FALSE])
  rownames(cells) <- NULL
  cells
}

#' rows to run. n up to MAX_N. w/ CELLS=ends only the first and last per spec and fold count
cells_to_run <- function(grid) {
  cells <- if (is.null(grid$n)) seq_len(nrow(grid)) else which(grid$n <= as.numeric(Sys.getenv("MAX_N", Inf)))
  if (Sys.getenv("CELLS") != "ends") return(cells)
  by <- intersect(c("spec", "folds"), names(grid))
  version <- if (length(by)) do.call(paste, grid[by]) else rep("", nrow(grid))
  sort(unname(unlist(lapply(split(cells, version[cells]), function(v) unique(c(head(v, 1), tail(v, 1)))))))
}

# study 2 extra step. a pilot lambda grid per cell on its own seed (Mersenne-Twister 1000 + block)

#' every cell's pilot lambda grid. also written to <grid>_lamgrid.csv (design cols block k lambda)
pilot_grids <- function(grid, generate, settings, out_file, n_lambda = 30) {
  block  <- stream_blocks(grid)
  design <- grid[setdiff(names(grid), "block")]
  grids  <- lapply(seq_len(nrow(grid)), function(i) {
    cell <- as.list(design[i, , drop = FALSE])
    set.seed(1000 + block[i], kind = "Mersenne-Twister")
    s <- settings(cell)
    pilot_lambda_grid(generate(cell), n_lambda, s$alpha, s$max_imbalance, s$min_ratio)
  })
  rows <- rep(seq_len(nrow(grid)), lengths(grids))
  data.table::fwrite(data.table::data.table(design[rows, , drop = FALSE], block = block[rows],
                                            k = sequence(lengths(grids)), lambda = unlist(grids)),
                     sub("\\.csv\\.gz$", "_lamgrid.csv", out_file))
  grids
}

#' one study 2 grid. pilot grids then estimate_study2() per cell. 10 folds unless a folds column
run_study2_grid <- function(generate, grid, replications, base_seed, out_file,
                            settings = cell_path_settings, workers = 1,
                            stream_block = replications, cells = seq_len(nrow(grid))) {
  if (file.exists(out_file)) {
    message(basename(out_file), " exists, skipped")
    return(invisible(NULL))
  }
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  design <- grid[setdiff(names(grid), "block")]
  key    <- function(cell) paste(unlist(cell), collapse = " ")
  grids  <- setNames(pilot_grids(grid, generate, settings, out_file),
                     vapply(seq_len(nrow(design)), function(i) key(as.list(design[i, , drop = FALSE])), ""))
  with_grid <- function(cell) c(generate(cell), list(lambda_grid = grids[[key(cell)]], cell = cell))
  estimate  <- function(d) {
    s <- settings(d$cell)
    estimate_study2(d, d$lambda_grid, nfolds = if (is.null(d$cell$folds)) 10 else d$cell$folds,
                    alpha = s$alpha, max_imbalance = s$max_imbalance, min_ratio = s$min_ratio)
  }
  run_grid(with_grid, grid, replications, base_seed, out_file, estimate, workers, stream_block, cells)
}

#' a cell's penalty and path end. lasso at 1e-4 unless the grid has alpha and min_ratio columns
cell_path_settings <- function(cell) {
  alpha <- if (is.null(cell$alpha)) 1 else cell$alpha
  list(alpha = alpha, max_imbalance = 1e-4,
       min_ratio = if (!is.null(cell$min_ratio)) cell$min_ratio else if (alpha == 0) 1e-4 else 1e-2)
}
