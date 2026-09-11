# Build cleaned summary tables for the CV exploration ----
library(dplyr)
dir <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/Exploring CV"
stopifnot(dir.exists(dir))
sels <- c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")

# families to keep (>= 200 reps, identified design); alpha for those whose
# cell has no alpha column
keep_fam <- c("tune", "tune2", "tune4", "tunea3", "tunea4", "tunen2", "tunep2",
              "snr_enet", "n_enet", "alpha_bad", "spread", "dima", "dimhi",
              "dimsnr", "snr_good", "snr_mb", "ov_s1", "ov_1k")
alpha_of <- c(tune = 1, tune2 = 1, tune4 = 1, tunen2 = 1, tunep2 = 1,
              snr_enet = 0.5, n_enet = 0.5, snr_good = 1, snr_mb = 1,
              ov_s1 = 1, ov_1k = 1)
fam_label <- c(tune = "SNR x overlap, diluted confounding (lasso)",
               tune2 = "SNR x overlap, fixed confounding (lasso)",
               tune4 = "SNR x overlap (lasso)",
               snr_enet = "SNR x overlap (elastic net)",
               tunea3 = "overlap x density, ridge/EN",
               tunea4 = "overlap, ridge/EN, floor 1e-4",
               alpha_bad = "alpha sweep, bad overlap, SNR",
               tunen2 = "n x overlap x SNR (lasso)",
               n_enet = "n x overlap x SNR (EN, bad/awful)",
               tunep2 = "p x overlap x density (lasso, 50 reps)",
               spread = "confounder spread x overlap x penalty",
               dima = "p x alpha x overlap",
               dimhi = "p in {1000, 2000} x alpha x overlap",
               dimsnr = "SNR at high p x alpha x overlap",
               snr_good = "noise axis, good overlap (lasso)",
               snr_mb = "noise axis, moderate/bad, sigma_y 4 and 7 (lasso)",
               ov_s1 = "overlap axis, sigma_y = 1 (lasso, maxit 1e4)",
               ov_1k = "overlap axis, sigma_y = 1 (lasso, maxit 1e5, 1000 reps)")


fams <- sub("_?\\d+_r\\d+\\.rds$", "", list.files(dir, "_r\\d+\\.rds$"))
table(fams)[!names(table(fams)) %in% keep_fam]

summarise_file <- function(f) {
  x <- readRDS(f)
  ok <- vapply(x$res, \(r) is.null(r$err), logical(1))
  res <- x$res[ok]
  reached <- if (is.null(res[[1]]$lam_end)) rep(1, length(x$lam)) else
    colMeans(outer(vapply(res, `[[`, numeric(1), "lam_end"), x$lam, "<="))
  keep <- reached >= 0.95
  cell <- x$cell
  sig <- if (!is.null(x$sigmas)) x$sigmas else
    if (!is.null(cell$sigma_y)) cell$sigma_y else 1
  cell$sigma_y <- NULL
  by_sigma <- lapply(seq_along(sig), function(s) {
    if (is.null(x$sigmas)) {
      ep <- do.call(rbind, lapply(res, `[[`, "est_path"))
      es <- do.call(rbind, lapply(res, `[[`, "est_sel"))
    } else {
      ep <- do.call(rbind, lapply(res, \(r) r$est_path[, s]))
      es <- do.call(rbind, lapply(res, \(r) r$est_sel[s, ]))
    }
    rp <- sqrt(colMeans(ep^2))
    rp[!keep] <- NA
    j_min <- which.min(rp)
    j_floor <- max(which(keep))
    rs <- sqrt(colMeans(es^2))[sels]
    # Mean selected lambdas identify the same cases as the existing pick tables.
    # Shared-noise files may store one lambda vector or one row per sigma.
    lm <- colMeans(do.call(rbind, lapply(res, function(r) {
      v <- r$lam_sel
      if (is.null(v)) return(c(cv.bloss = NA_real_, boot.inf = NA_real_))
      if (!is.null(dim(v))) v <- v[s, ]
      v[c("cv.bloss", "boot.inf")]
    })))
    best <- sels[which.min(rs)]
    z_vs_floor <- function(k) {                  # negative = selector better
      d <- es[, k]^2 - ep[, j_floor]^2
      mean(d) / (sd(d) / sqrt(length(d)))
    }
    data.frame(file = basename(f), cell, sigma_y = sig[s], n_ok = sum(ok),
               lam_min = x$lam[j_min], lam_floor = x$lam[j_floor],
               rmse_min = rp[j_min], rmse_floor = rp[j_floor],
               gain_floor = rp[j_floor] / rp[j_min] - 1,
               as.list(setNames(rs / rp[j_min], paste0(sels, "_vs_min"))),
               lam_bloss = unname(lm["cv.bloss"]),
               lam_bootinf = unname(lm["boot.inf"]),
               # True target is zero; es contains each replication's selected estimate.
               RMSE_cv = unname(rs["cv.bloss"]),
               RMSE_boot = unname(rs["boot.inf"]),
               CV_excess_pct = 100 * (unname(rs["cv.bloss"]) / rp[j_min] - 1),
               Boot_excess_pct = 100 * (unname(rs["boot.inf"]) / rp[j_min] - 1),
               best = best, rmse_best = rs[best],
               best_vs_floor = rs[best] / rp[j_floor] - 1,
               z_best = z_vs_floor(best), z_bootinf = z_vs_floor("boot.inf"),
               z_cvbloss = z_vs_floor("cv.bloss"))
  })
  do.call(rbind, by_sigma)
}

files <- list.files(dir, "_r\\d+\\.rds$", full.names = TRUE)
tab <- bind_rows(lapply(files, summarise_file)) |>
  mutate(family = sub("_?\\d+_r\\d+\\.rds$", "", file)) |>
  filter(family %in% keep_fam, n_ok >= 200) |>
  group_by(family, across(any_of(c("overlap", "c_prop", "sigma_y", "n", "p",
                                   "s", "s_y", "dense", "alpha")))) |>
  slice_max(n_ok, n = 1, with_ties = FALSE) |>   # largest rep file per cell
  ungroup() |>
  mutate(alpha = coalesce(alpha, unname(alpha_of[family])),
         design = unname(fam_label[family]),
         r2 = 2.23 / (2.23 + sigma_y^2),
         mcse_rel = 1 / sqrt(2 * n_ok)) |>
  relocate(family, design, alpha)
readr::write_csv(tab, file.path(dir, "cv_summary.csv"))

# win counts per family x alpha (z beyond 2 = win; interior = floor > 5% above min)
counts <- tab |>
  group_by(family, design, alpha) |>
  summarise(cells = n(), interior = sum(gain_floor > 0.05),
            bootinf_wins = sum(z_bootinf < -2),
            cvbloss_wins = sum(z_cvbloss < -2),
            any_sel_win = sum(z_best < -2), floor_wins = sum(z_best > 2),
            .groups = "drop")
readr::write_csv(counts, file.path(dir, "cv_family_counts.csv"))

# best EN/ridge selector against the lasso floor for the same design
vs_floor <- function(d, by) {
  d |>
    group_by(across(all_of(by))) |>
    filter(any(alpha == 1), any(alpha < 1)) |>
    summarise(lasso_floor = rmse_floor[alpha == 1][1],
              best_sel = min(rmse_best[alpha < 1]),
              best = best[alpha < 1][which.min(rmse_best[alpha < 1])],
              best_alpha = alpha[alpha < 1][which.min(rmse_best[alpha < 1])],
              ratio = best_sel / lasso_floor, .groups = "drop")
}
vs_lasso <- bind_rows(
  vs_floor(filter(tab, family %in% c("tune4", "snr_enet")),
           c("overlap", "sigma_y")) |> mutate(design = "n = 1000, p = 100"),
  vs_floor(filter(tab, family %in% c("tunen2", "n_enet")),
           c("overlap", "n", "sigma_y")) |> mutate(design = "n sweep, p = 100"),
  vs_floor(filter(tab, family == "alpha_bad"), c("sigma_y")) |>
    mutate(design = "alpha sweep, bad overlap", overlap = "bad"),
  vs_floor(filter(tab, family == "spread"), c("overlap", "s")) |>
    mutate(design = "confounder spread", sigma_y = 1),
  vs_floor(filter(tab, family %in% c("dimhi", "dimsnr")),
           c("overlap", "p", "sigma_y")) |> mutate(design = "high p")) |>
  relocate(design, overlap, sigma_y, n, p, s)
readr::write_csv(vs_lasso, file.path(dir, "cv_vs_lasso_floor.csv"))
print(vs_lasso, n = Inf)

# Rebuild RMSE tables from simulation cases in the current summary.
# Use a single family for each table. ov_1k is the 1000-rep overlap sweep;
# use ov_s1 only when ov_1k is absent from the saved summary.
noise_family <- "snr_good"
overlap_family <- if ("ov_1k" %in% tab$family) "ov_1k" else "ov_s1"

make_rmse_table <- function(d, axis) {
  if (nrow(d) == 0L) {
    stop("No summary rows for ", axis, ". Available families: ",
         paste(sort(unique(tab$family)), collapse = ", "))
  }
  if (anyNA(d[[axis]]) || anyDuplicated(d[[axis]])) {
    print(d |> dplyr::select(dplyr::any_of(
      c("file", "family", axis, "n", "p", "overlap", "c_prop", "sigma_y", "alpha"))))
    stop(axis, " must identify one simulation case per row; ",
         "filter the design shown above before exporting.")
  }
  d |>
    dplyr::arrange(.data[[axis]]) |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(axis)),
      lam_opt = lam_min,
      dplyr::across(dplyr::any_of(c("lam_bloss", "lam_bootinf"))),
      rmse_min,
      RMSE_cv = rmse_min * .data[["cv.bloss_vs_min"]],
      RMSE_boot = rmse_min * .data[["boot.inf_vs_min"]],
      CV_excess_pct = 100 * (.data[["cv.bloss_vs_min"]] - 1),
      Boot_excess_pct = 100 * (.data[["boot.inf_vs_min"]] - 1),
      source_file = file, n_ok
    )
}

noise_rmse <- make_rmse_table(
  dplyr::filter(tab, family == noise_family, alpha == 1), "sigma_y")
overlap_rmse <- make_rmse_table(
  dplyr::filter(tab, family == overlap_family, alpha == 1, sigma_y == 1), "c_prop")

# Compute both tables successfully before writing either output.
readr::write_csv(noise_rmse, file.path(dir, "cv_rmse_noise.csv"))
readr::write_csv(overlap_rmse, file.path(dir, "cv_rmse_overlap.csv"))
message("Exported noise family: ", noise_family,
        "; overlap family: ", overlap_family)





