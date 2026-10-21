# performance measures (Methods 3.5) from results files. tau = 0 so every error is the estimate
# a cell = design level x outcome (x spec where pooled). MCSEs per Morris et al. 2019

library(dplyr)
library(tidyr)

z_95 <- qnorm(0.975)

#' one results file in long form (design cols level sim estimator outcome error)
read_results <- function(file, estimators = NULL) {
  wide <- data.table::fread(file, colClasses = list(character = "err"))
  if (!is.null(estimators)) wide <- wide[wide$estimator %in% estimators, ]
  if ("batch" %in% names(wide)) {                 # the n = 1 million run. one cell in 10 batches
    wide$sim <- (wide$batch - 1L) * max(wide$sim) + wide$sim
    wide$batch <- NULL
  }
  outcome_cols <- intersect(c("linear", "quadratic", "exponential", "quad1", "exp"), names(wide))
  design_cols  <- setdiff(names(wide), c("sim", "estimator", "err", outcome_cols))
  as_text <- function(v) if (is.numeric(v))        # 100000 and 1e+05 give the same level
    format(v, scientific = FALSE, drop0trailing = TRUE, trim = TRUE, digits = 15) else as.character(v)
  as_tibble(wide) |>
    mutate(level = do.call(paste, c(lapply(across(all_of(design_cols)), as_text), sep = "_"))) |>
    select(all_of(design_cols), level, sim, estimator, all_of(outcome_cols)) |>
    pivot_longer(all_of(outcome_cols), names_to = "outcome", values_to = "error") |>
    mutate(error = suppressWarnings(as.numeric(error)),
           outcome = recode(outcome, quad1 = "quadratic", exp = "exponential"))
}

# per rep influence values for the delta method MCSEs of R_c s Rbar and sbar
influence_log_mse <- function(e) e^2 / mean(e^2) - 1
influence_share   <- function(e) {
  bias <- mean(e); mse <- mean(e^2)
  (2 * bias * (e - bias) - bias^2 / mse * (e^2 - mse)) / mse
}

#' one estimator in one cell. coverage = share w/ |error| <= 1.96 EmpSE. s = Bias^2 / RMSE^2
cell_measures <- function(error) {
  e <- error[is.finite(error)]
  m <- length(e)
  if (m < 2) return(tibble(reps = length(error), missing = length(error) - m, m = m))
  rmse  <- sqrt(mean(e^2))
  empse <- sd(e)
  cover <- mean(abs(e) <= z_95 * empse)
  tibble(reps = length(error), missing = length(error) - m, m = m,
         bias = mean(e),       mcse_bias     = empse / sqrt(m),
         empse = empse,        mcse_empse    = empse / sqrt(2 * (m - 1)),
         rmse = rmse,          mcse_rmse     = sd(e^2) / (2 * rmse * sqrt(m)),
         coverage = cover,     mcse_coverage = sqrt(cover * (1 - cover) / m),
         share = mean(e)^2 / mean(e^2), mcse_share = sd(influence_share(e)) / sqrt(m))
}

#' every estimator in every cell. dropped = cell missing over 1% of its estimates
performance <- function(results, max_missing = 0.01) {
  results |>
    summarise(cell_measures(error), .by = all_of(setdiff(names(results), c("sim", "error")))) |>
    mutate(dropped = missing / reps > max_missing)
}

#' A and B errors on the same draws w/ influence values. a cell is dropped if either misses over 1%
#' below that the reps missing either estimate are dropped from both
pair_draws <- function(results, A, B, max_missing = 0.01) {
  wide <- results |>
    filter(estimator %in% c(A, B)) |>
    select(spec, level, sim, outcome, estimator, error) |>
    pivot_wider(names_from = estimator, values_from = error) |>
    rename(e_A = all_of(A), e_B = all_of(B))
  cells <- wide |>
    summarise(miss_A = mean(!is.finite(e_A)), miss_B = mean(!is.finite(e_B)), .by = c(spec, level, outcome)) |>
    mutate(kept = miss_A <= max_missing & miss_B <= max_missing)
  draws <- wide |>
    semi_join(filter(cells, kept), by = c("spec", "level", "outcome")) |>
    filter(is.finite(e_A), is.finite(e_B)) |>
    mutate(influence_log_R = (influence_log_mse(e_A) - influence_log_mse(e_B)) / 2,
           influence_s_A = influence_share(e_A), influence_s_B = influence_share(e_B),
           .by = c(spec, level, outcome))
  list(draws = draws, cells = cells)
}

#' R_c and both shares of bias per kept cell. won and lost = 3 MCSEs from zero on the log scale
pair_cells <- function(draws) {
  draws |>
    summarise(reps = n(), rmse_A = sqrt(mean(e_A^2)), rmse_B = sqrt(mean(e_B^2)),
              R_c = rmse_A / rmse_B, log_R_c_mcse = sd(influence_log_R) / sqrt(reps),
              s_A = mean(e_A)^2 / mean(e_A^2), s_A_mcse = sd(influence_s_A) / sqrt(reps),
              s_B = mean(e_B)^2 / mean(e_B^2), s_B_mcse = sd(influence_s_B) / sqrt(reps),
              .by = c(spec, level, outcome)) |>
    mutate(won = log(R_c) + 3 * log_R_c_mcse < 0, lost = log(R_c) - 3 * log_R_c_mcse > 0)
}

#' Rbar = 100 (1 - exp(mean log R_c)) and sbar over the given cells w/ the shared draw MCSE
#' a draw's influence values are summed over the cells it feeds before the variance. per level
average_pair <- function(draws, cells = pair_cells(draws)) {
  if (nrow(cells) == 0) return(tibble(n_cells = 0L))
  v <- draws |>
    semi_join(cells, by = c("spec", "level", "outcome")) |>
    summarise(across(c(influence_log_R, influence_s_A, influence_s_B), sum), .by = c(level, sim)) |>
    summarise(across(-sim, \(g) var(g) / n()), .by = level) |>
    summarise(across(-level, sum))
  mcse <- sqrt(v) / nrow(cells)
  mean_log_R <- mean(log(cells$R_c))
  tibble(n_cells = nrow(cells),
         R_bar = 100 * (1 - exp(mean_log_R)), R_bar_mcse = 100 * exp(mean_log_R) * mcse$influence_log_R,
         R_min = min(cells$R_c), R_max = max(cells$R_c), won = sum(cells$won), lost = sum(cells$lost),
         s_bar_A = mean(cells$s_A), s_bar_A_mcse = mcse$influence_s_A,
         s_bar_B = mean(cells$s_B), s_bar_B_mcse = mcse$influence_s_B)
}

#' Rbar and sbar of A vs B. pooled and by outcome. per spec and over both
relative_rmse <- function(results, A, B, max_missing = 0.01) {
  pd    <- pair_draws(results, A, B, max_missing)
  cells <- pair_cells(pd$draws)
  blocks <- unique(cells$spec)
  if (length(blocks) > 1) blocks <- c(blocks, "both")
  lapply(blocks, function(block) {
    in_block <- if (block == "both") cells else filter(cells, spec == block)
    by_outcome <- lapply(group_split(in_block, outcome), \(k) mutate(average_pair(pd$draws, k), outcome = k$outcome[1]))
    bind_rows(mutate(average_pair(pd$draws, in_block), outcome = "pooled"), by_outcome) |>
      mutate(block = block, .before = 1)
  }) |>
    bind_rows() |>
    mutate(A = A, B = B, .before = 1)
}

#' study 2 rules per cell. interior = path min more than 5% below the endpoint RMSE
#' win = a selector RMSE below the endpoint by over 3 combined MCSEs. lead = lowest RMSE selector
study2_rules <- function(perf, keys = c("spec", "level", "outcome", "strength", "est"), gain = 0.05, z = 3) {
  keys <- intersect(keys, names(perf))
  end <- perf |> filter(arm == "fixed") |> select(all_of(keys), reps, end = rmse, end_mcse = mcse_rmse)
  path_min <- perf |> filter(arm == "path") |> summarise(path_min = min(rmse), .by = all_of(keys))
  selectors <- perf |>
    filter(arm %in% c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")) |>
    inner_join(end, by = keys) |>
    summarise(win = any(rmse < end - z * sqrt(mcse_rmse^2 + end_mcse^2)), lead = arm[which.min(rmse)],
              .by = all_of(keys))
  end |>
    inner_join(path_min, by = keys) |>
    inner_join(selectors, by = keys) |>
    mutate(interior = path_min < (1 - gain) * end)
}
