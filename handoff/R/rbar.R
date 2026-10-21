# Study 1 performance measures (thesis Performance measures) and averages over cells.
# only place they get computed. build.R stores them verify.R rechecks separately thesis.R reads

library(dplyr)
library(tidyr)
library(purrr)

# influence values (IF). smooth fns of means so delta method MCSE = sd(IF_r) / sqrt(m)
# kept per rep since averages over cells that share draws need them

#' IF of log MSE. finite errors only
influence_log_mse <- function(e) e^2 / mean(e^2) - 1

#' IF of the share of bias s = Bias^2 / MSE
influence_share <- function(e) {
  bias <- mean(e)
  mse <- mean(e^2)
  (2 * bias * (e - bias) - bias^2 / mse * (e^2 - mse)) / mse
}

# one estimator in one cell. tau = 0 in all Study 1 cells so error = stored estimate

#' bias empse rmse and share of bias of one estimator in one cell w/ MCSEs. NA reps left out
#' n_extreme = count of |e_r| > 10. heavy tails make the delta MCSE optimistic
arm_measures <- function(error) {
  e <- error[is.finite(error)]
  m <- length(e)
  counts <- tibble(reps = length(error), missing = length(error) - m, m = m,
                   n_extreme = sum(abs(e) > 10))
  if (m < 2) {
    na <- NA_real_
    return(mutate(counts, bias = na, empse = na, rmse = na, mcse_bias = na,
                  mcse_empse = na, mcse_rmse = na, share = na,
                  mcse_share = na))
  }
  rmse <- sqrt(mean(e^2))
  mutate(counts,
         bias = mean(e),
         empse = sd(e),
         rmse = rmse,
         mcse_bias = sd(e) / sqrt(m),
         mcse_empse = sd(e) / sqrt(2 * (m - 1)),
         mcse_rmse = sd(e^2) / (2 * rmse * sqrt(m)),
         share = mean(e)^2 / mean(e^2),
         mcse_share = sd(influence_share(e)) / sqrt(m))
}

#' performance of each estimator in each cell. dropped = over miss_max (1%) of estimates missing
#' tables show a dash there and the pair averages leave it out
cell_performance <- function(est, miss_max = 0.01) {
  est |>
    summarise(arm_measures(error), .by = c(spec, level, estimator, outcome)) |>
    mutate(dropped = missing / reps > miss_max)
}

# estimator A vs comparator B in one cell

#' errors of A and B on the same draws + IFs. pairing makes R_c far more precise than 2 indep RMSEs
#' cell out if either misses > miss_max else pairwise deletion. 1% is our rule
pair_draws <- function(est, A, B, miss_max = 0.01) {
  wide <- est |>
    filter(estimator %in% c(A, B)) |>
    select(spec, level, sim, outcome, estimator, error) |>
    pivot_wider(names_from = estimator, values_from = error) |>
    rename(e_A = all_of(A), e_B = all_of(B))
  cells <- wide |>
    summarise(miss_A = mean(!is.finite(e_A)), miss_B = mean(!is.finite(e_B)),
              .by = c(spec, level, outcome)) |>
    mutate(kept = miss_A <= miss_max & miss_B <= miss_max)
  draws <- wide |>
    semi_join(filter(cells, kept), by = c("spec", "level", "outcome")) |>
    filter(is.finite(e_A), is.finite(e_B)) |>
    mutate(
      # log R_c = (log MSE_A - log MSE_B) / 2 so its IF is half the difference
      influence_log_R = (influence_log_mse(e_A) - influence_log_mse(e_B)) / 2,
      influence_s_A = influence_share(e_A),
      influence_s_B = influence_share(e_B),
      .by = c(spec, level, outcome))
  list(draws = draws, cells = cells)
}

#' R_c both shares of bias and the paired RMSE diff in each kept cell w/ MCSEs
#' won / lost = log R_c over 3 MCSE below / above 0. no multiplicity adj
pair_cell_measures <- function(draws) {
  draws |>
    summarise(
      reps = n(),
      rmse_A = sqrt(mean(e_A^2)),
      rmse_B = sqrt(mean(e_B^2)),
      R_c = rmse_A / rmse_B,
      log_R_c_mcse = sd(influence_log_R) / sqrt(reps),
      # paired RMSE diff for tbl-q1-ols. IF of an RMSE is e_r^2 / (2 RMSE)
      rmse_diff = rmse_A - rmse_B,
      rmse_diff_mcse = sd(e_A^2 / (2 * rmse_A) - e_B^2 / (2 * rmse_B)) /
        sqrt(reps),
      s_A = mean(e_A)^2 / mean(e_A^2),
      s_A_mcse = sd(influence_s_A) / sqrt(reps),
      s_B = mean(e_B)^2 / mean(e_B^2),
      s_B_mcse = sd(influence_s_B) / sqrt(reps),
      .by = c(spec, level, outcome)) |>
    mutate(won = log(R_c) + 3 * log_R_c_mcse < 0,
           lost = log(R_c) - 3 * log_R_c_mcse > 0)
}

# averages over Theta. equal weight means. Rbar on log scale so a few big cells can't dominate
# cells share draws so IFs summed per draw (g_r) first. indep would understate MCSE 1.3-1.8x

#' IFs of a pair by design level. one draw x cell matrix per measure. 0 where no paired est
level_matrices <- function(draws) {
  draws <- mutate(draws, cell = paste(spec, outcome, sep = "|"))
  map(split(draws, draws$level), \(d) {
    sims <- sort(unique(d$sim))
    cells <- unique(d$cell)
    at <- cbind(match(d$sim, sims), match(d$cell, cells))
    as_matrix <- \(values) {
      M <- matrix(0, length(sims), length(cells), dimnames = list(NULL, cells))
      M[at] <- values
      M
    }
    list(present = as_matrix(1),
         influence_log_R = as_matrix(d$influence_log_R),
         influence_s_A = as_matrix(d$influence_s_A),
         influence_s_B = as_matrix(d$influence_s_B))
  })
}

#' MCSEs of mean log R_c sbar_A sbar_B over Theta. levels have own seed streams so variances add
mean_mcse_over_cells <- function(by_level, cells) {
  v <- c(influence_log_R = 0, influence_s_A = 0, influence_s_B = 0)
  for (lv in unique(cells$level)) {
    at <- cells$level == lv
    cols <- paste(cells$spec[at], cells$outcome[at], sep = "|")
    level <- by_level[[lv]]
    draws <- rowSums(level$present[, cols, drop = FALSE]) > 0  # m_level
    for (f in names(v)) {
      g <- rowSums(level[[f]][draws, cols, drop = FALSE])     # g_r
      v[f] <- v[f] + var(g) / sum(draws)
    }
  }
  sqrt(v) / nrow(cells)
}

#' Rbar and sbar of one pair over Theta w/ MCSEs + range of R_c and won / lost counts
average_over_cells <- function(by_level, cells) {
  k <- nrow(cells)
  if (k == 0) {
    na <- NA_real_
    return(tibble(n_cells = 0L, mean_log_R = na, mean_log_R_mcse = na,
                  R_bar = na, R_bar_mcse = na, R_min = na, R_max = na,
                  won = 0L, lost = 0L, s_bar_A = na, s_bar_A_mcse = na,
                  s_bar_B = na, s_bar_B_mcse = na))
  }
  se <- as.list(mean_mcse_over_cells(by_level, cells))
  mean_log_R <- mean(log(cells$R_c))
  tibble(n_cells = k,
         mean_log_R = mean_log_R,
         mean_log_R_mcse = se$influence_log_R,
         R_bar = 100 * (1 - exp(mean_log_R)),
         R_bar_mcse = 100 * exp(mean_log_R) * se$influence_log_R,
         R_min = min(cells$R_c),
         R_max = max(cells$R_c),
         won = sum(cells$won),
         lost = sum(cells$lost),
         s_bar_A = mean(cells$s_A),
         s_bar_A_mcse = se$influence_s_A,
         s_bar_B = mean(cells$s_B),
         s_bar_B_mcse = se$influence_s_B)
}

# appendix measures. MAE kurtosis kurtosis-adjusted EmpSE MCSE and approx bias vs OLS

#' bootstrap SE of the median. each resample median drawn from Beta order stats of uniforms
#' same bootstrap distribution but 2 draws per resample not m
median_se_bootstrap <- function(a, B = 1000) {
  m <- length(a)
  if (m < 2 || B < 2) return(NA_real_)
  a <- sort(a)
  k <- floor((m + 1) / 2)                   # lower middle order stat
  at <- \(u) a[pmin(floor(m * u) + 1, m)]
  u1 <- rbeta(B, k, m - k + 1)
  med <- if (m %% 2 == 1) at(u1) else {
    u2 <- u1 + (1 - u1) * rbeta(B, 1, m - k)
    (at(u1) + at(u2)) / 2
  }
  sd(med)
}

#' arm_measures() + kurtosis and MAE (median |e_r| as Kang and Schafer 2007) w/ bootstrap MCSE
#' mcse_empse_k allows for kurtosis as Morris et al 2019 Table 6 assumes normal estimates
arm_measures_full <- function(error, B = 1000, extreme_at = 10) {
  base <- arm_measures(error)
  e <- error[is.finite(error)]
  base$n_extreme <- sum(abs(e) > extreme_at)
  if (length(e) < 2) {
    na <- NA_real_
    return(mutate(base, kurt = na, mcse_empse_k = na, mae = na,
                  mcse_mae = na))
  }
  centred <- e - mean(e)
  kurt <- mean(centred^4) / mean(centred^2)^2 - 3
  mutate(base,
         kurt = kurt,
         mcse_empse_k = mcse_empse * sqrt(1 + kurt / 2),
         mae = median(abs(e)),
         mcse_mae = median_se_bootstrap(abs(e), B))
}

#' bias of A and B and approx bias (paired mean of e_A - e_B) in each kept cell w/ MCSEs
#' w/ A a balancing arm and B = OLS it shows how far the weights are from the regression
approx_cell_measures <- function(draws) {
  draws |>
    summarise(reps = n(),
              bias_A = mean(e_A), mcse_bias_A = sd(e_A) / sqrt(reps),
              bias_B = mean(e_B), mcse_bias_B = sd(e_B) / sqrt(reps),
              bias_diff = mean(e_A - e_B),
              mcse_bias_diff = sd(e_A - e_B) / sqrt(reps),
              .by = c(spec, level, outcome))
}

#' MCSE of a mean over cells of one IF column. sum per draw then var within level then add levels
mean_mcse_of <- function(draws, influence, k) {
  draws |>
    summarise(g = sum(.data[[influence]]), .by = c(level, sim)) |>
    summarise(v = var(g) / n(), .by = level) |>
    pull(v) |>
    sum() |>
    sqrt() / k
}

#' approx bias averaged over cells (equal weight and signed) pooled then per outcome model
#' bias of A and B beside it. IF of a cell mean is the centred error
approx_average <- function(draws) {
  draws <- mutate(draws,
                  influence_A = e_A - mean(e_A), influence_B = e_B - mean(e_B),
                  influence_diff = (e_A - e_B) - mean(e_A - e_B),
                  .by = c(spec, level, outcome))
  cells <- summarise(draws, bias_A = mean(e_A), bias_B = mean(e_B),
                     bias_diff = mean(e_A - e_B), .by = c(spec, level, outcome))
  one_set <- \(d, cl) {
    k <- nrow(cl)
    tibble(n_cells = k,
           bias_A = mean(cl$bias_A),
           mcse_bias_A = mean_mcse_of(d, "influence_A", k),
           bias_B = mean(cl$bias_B),
           mcse_bias_B = mean_mcse_of(d, "influence_B", k),
           bias_diff = mean(cl$bias_diff),
           mcse_bias_diff = mean_mcse_of(d, "influence_diff", k))
  }
  models <- intersect(c("linear", "quad1", "exp"), cells$outcome)
  bind_rows(
    mutate(one_set(draws, cells), outcome = "pooled"),
    map(models, \(o) {
      mutate(one_set(filter(draws, outcome == o), filter(cells, outcome == o)),
             outcome = o)
    }) |>
      list_rbind()
  )
}
