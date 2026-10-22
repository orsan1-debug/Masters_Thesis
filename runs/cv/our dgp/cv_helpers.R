# shared steps of the cv run files in this folder
# run file sources R/dgp.R estimators_ipw.R estimators_cv_dgp2.R simulate.R path_floor.R first

#' outcome strengths of the cv grids. halving from 1 then 0
cv_strengths <- c(1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0)

#' one dgp2() draw for the cv runs. same setup as runs/aipw/overlap_axis.R
#' misspec = analyst sees Kang-Schafer transforms of the 4 ps covs. Y and m still from X_true
dgp2_cv <- function(overlap, n = 1000, p = 50, misspec = FALSE, strength = 1)
  dgp2(n = n, p = p, s = 4, signs = "pos",
       outcome = c("linear", "quad1", "exp"), misspec = misspec,
       covcor = "iid", overlap = overlap, strength = strength)

#' every strength on the same draw. Y = s * m + eps per strength (eps = Y - m as tau = 0)
#' cols eg linear_s0.125. all 7 strengths share the weights folds and boot half samples
share_strengths <- function(d, strengths = cv_strengths) {
  eps <- d$Y - d$m
  d$Y <- do.call(cbind, lapply(strengths, function(s)
    `colnames<-`(s * d$m + eps, paste0(colnames(d$m), "_s", s))))
  d
}

#' lambda grid per cell from a pilot draw on MT seed 1000 + cell so a resumed run rebuilds it
#' written next to the results as <stem>_lamgrid.csv
pilot_lam_grids <- function(grid, dgp_gen, out_file) {
  lam <- lapply(grid$cell, function(i) {
    set.seed(1000 + i, kind = "Mersenne-Twister")
    f <- cell_floor(grid[i, ])
    make_lam_grid(dgp_gen(grid[i, ]), alpha = f$alpha,
                  max_imbalance = f$max_imbalance, min_ratio = f$min_ratio)
  })
  data.table::fwrite(data.frame(cell = rep(grid$cell, lengths(lam)),
                                k = sequence(lengths(lam)),
                                lambda = unlist(lam)),
                     sub("\\.csv\\.gz$", "_lamgrid.csv", out_file))
  lam
}

#' draw plus lambda grid and cell. factory w/ own names as future rebinds variables by name
#' a closure calling dgp_gen inside simulate_grid() called itself till R stopped it
cv_dgp_gen <- function(cv_draw, cv_lam) {
  force(cv_draw)
  force(cv_lam)
  function(cell) c(cv_draw(cell), list(lam_grid = cv_lam[[cell$cell]],
                                       cell = cell))
}

#' estimate_cv() at the cell's penalty and path floor
cv_estimator <- function(cv_nfolds) {
  force(cv_nfolds)
  function(d) {
    f <- cell_floor(d$cell)
    estimate_cv(d, nfolds = cv_nfolds, alpha = f$alpha,
                max_imbalance = f$max_imbalance, min_ratio = f$min_ratio)
  }
}

#' run one cv grid. skips if results file exists. on error reruns from done cells up to tries times
#' forked workers under Rscript on mac/linux. windows keeps multisession (no multicore)
run_cv_grid <- function(dgp_gen, grid, num_sim, base_seed, out_file,
                        nfolds = 10, tries = 5) {
  if (file.exists(out_file)) {
    message(basename(out_file), " exists: skipped")
    return(invisible(out_file))
  }
  stopifnot(identical(grid$cell, seq_len(nrow(grid))))
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  if (is.null(getOption("sim.plan")) && !interactive() &&
      .Platform$OS.type == "unix") {
    old <- options(sim.plan = future::multicore)
    on.exit(options(old), add = TRUE)
  }

  lam <- pilot_lam_grids(grid, dgp_gen, out_file)
  ck_dir <- sub("\\.csv\\.gz$", "_cells", out_file)
  for (attempt in seq_len(tries)) {
    message(sprintf("%s, try %d: %d/%d cells done, %d reps each",
                    basename(out_file), attempt, length(list.files(ck_dir)),
                    nrow(grid), num_sim))
    ok <- tryCatch({
      simulate_grid(cv_dgp_gen(dgp_gen, lam), grid, num_sim = num_sim,
                    base_seed = base_seed, out_file = out_file,
                    estimator = cv_estimator(nfolds))
      TRUE
    }, error = function(e) {
      message("stopped: ", conditionMessage(e))
      FALSE
    })
    if (ok) return(invisible(out_file))
  }
  stop(basename(out_file), ": stopped ", tries, " times, see the messages above")
}
