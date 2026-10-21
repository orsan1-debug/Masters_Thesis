# synthetic tests of rbar.R w/ known answers + MC check of the Rbar and sbar MCSEs w/ shared draws
# no results files needed. ~1 min. source this file to run

source(here::here("handoff", "R", "rbar.R"))

#' one fake study in long form. m draws per level feed 3 outcomes w/ correlated A and B errors
make_study <- function(m = 500, levels = 2, sd_A = 0.6, sd_B = 1, b_A = 0.05,
                       b_B = 0.3, rho = 0.7) {
  out <- list()
  for (lv in seq_len(levels)) {
    z <- rnorm(m)                                    # shared by the draw
    for (o in c("linear", "quad1", "exp")) {
      u <- rho * z + sqrt(1 - rho^2) * rnorm(m)
      v <- rho * z + sqrt(1 - rho^2) * rnorm(m)
      eB <- b_B + sd_B * u
      eA <- b_A + sd_A * (0.8 * u + 0.6 * v)         # corr(A B) = 0.8
      out[[length(out) + 1]] <- tibble(
        spec = "correct", level = paste0("L", lv), sim = rep(seq_len(m), 2),
        outcome = o, estimator = rep(c("A", "B"), each = m),
        error = c(eA, eB))
    }
  }
  bind_rows(out)
}

#' Rbar and sbar of A vs B over all cells of a study via rbar.R
study_average <- function(est) {
  pd <- pair_draws(est, "A", "B")
  cells <- pair_cell_measures(pd$draws)
  average_over_cells(level_matrices(pd$draws), cells)
}

ok <- function(cond, what) {
  cat(if (isTRUE(cond)) "PASS" else "FAIL", what, "\n")
  invisible(isTRUE(cond))
}
results <- c()

# known answer. e_A = e_B / 2 in every cell so R_c = 0.5 and Rbar = 50
set.seed(1)
est <- make_study(m = 200)
est$error[est$estimator == "A"] <- est$error[est$estimator == "B"] / 2
a <- study_average(est)
results["known"] <- ok(abs(a$R_bar - 50) < 1e-12 && a$R_bar_mcse < 1e-12,
                       "e_A = e_B / 2: Rbar = 50, MCSE = 0")

# swap A and B. exp(mean log R) inverts so (100 - Rbar AB) * (100 - Rbar BA) = 100^2
# the percentages are not opposites
set.seed(2)
est <- make_study(m = 300)
ab <- study_average(est)
ba <- study_average(mutate(est, estimator = recode(estimator, A = "B",
                                                   B = "A")))
results["swap"] <- ok(abs((100 - ab$R_bar) * (100 - ba$R_bar) - 1e4) < 1e-8,
                      sprintf("swap: Rbar(A,B) = %.1f, Rbar(B,A) = %.1f",
                              ab$R_bar, ba$R_bar))

# jensen. R_c = 0.5 and 1.5 give Rbar = 13.4 not 0
results["jensen"] <- ok(abs(100 * (1 - exp(mean(log(c(0.5, 1.5))))) -
                              13.39746) < 1e-4,
                        "R_c = 0.5, 1.5: Rbar = 13.4 while mean(1 - R_c) = 0")

# bias share s = bias^2 / MSE = bias^2 / (bias^2 + (m - 1) / m * EmpSE^2)
e <- est$error[est$estimator == "B" & est$level == "L1" &
                 est$outcome == "linear"]
m <- length(e)
results["share"] <- ok(abs(mean(e)^2 / mean(e^2) -
                             mean(e)^2 / (mean(e)^2 + (m - 1) / m * var(e))) <
                         1e-12, "s = bias^2 / (bias^2 + (m-1)/m EmpSE^2)")

# mcse calibration. sd of Rbar and sbar over 400 studies vs mean delta MCSE. ratio near 1
set.seed(3)
sims <- replicate(400, unlist(study_average(make_study(m = 400))[
  c("R_bar", "R_bar_mcse", "s_bar_A", "s_bar_A_mcse", "s_bar_B",
    "s_bar_B_mcse")]))
ratio <- c(R_bar = sd(sims["R_bar", ]) / mean(sims["R_bar_mcse", ]),
           s_bar_A = sd(sims["s_bar_A", ]) / mean(sims["s_bar_A_mcse", ]),
           s_bar_B = sd(sims["s_bar_B", ]) / mean(sims["s_bar_B_mcse", ]))
results["mcse"] <- ok(all(ratio > 0.85 & ratio < 1.15),
                      paste("SD over studies / mean MCSE:",
                            paste(sprintf("%s %.2f", names(ratio), ratio),
                                  collapse = ", ")))
# for comparison. treating cells as independent (summing variances) understates it
naive <- replicate(200, {
  pd <- pair_draws(make_study(m = 400), "A", "B")
  cells <- pair_cell_measures(pd$draws)
  c(R = 100 * (1 - exp(mean(log(cells$R_c)))),
    se = 100 * exp(mean(log(cells$R_c))) *
      sqrt(sum(cells$log_R_c_mcse^2)) / nrow(cells))
})
cat(sprintf("  (cells treated as independent: SD / mean MCSE = %.2f)\n",
            sd(naive["R", ]) / mean(naive["se", ])))

# missing ests. a cell w/ > 1% missing is dropped from the pair
# below that draws missing either est are dropped from both arms
set.seed(4)
est <- make_study(m = 400)
idx <- which(est$estimator == "A" & est$level == "L1" & est$outcome == "exp")
est$error[idx[1:5]] <- NA       # 1.25% missing so cell dropped
idx2 <- which(est$estimator == "B" & est$level == "L2" &
                est$outcome == "exp")
est$error[idx2[1:2]] <- NA      # 0.5% missing so 2 draws dropped from both
pd <- pair_draws(est, "A", "B")
results["missing"] <- ok(
  sum(!pd$cells$kept) == 1 &&
    nrow(filter(pd$draws, level == "L2", outcome == "exp")) == 398,
  "1.25% missing drops the cell; 0.5% drops 2 draws from both arms")

# stop instead of print on a fail. the rebuild line sources the build right after
# this file and must not run if a test failed
if (!all(results)) {
  stop("Some tests FAILED: ", paste(names(results)[!results], collapse = ", "))
}
cat("All tests passed.\n")
