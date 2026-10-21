# path floor of the dgp2 cv runs in one place so no run relies on a function default
# floor = last lambda reached. estimate_cv() keeps it as the fixed arm and path row 30 repeats it

#' lasso floor. path ends at lambda = 1e-4 which bounds the largest |SMD| (balnet max.imbalance)
floor_max_imbalance <- 1e-4

#' lambda.min.ratio floor for alpha < 1. ridge 1e-4 since balnet refuses max.imbalance there
#' enet 1e-2 (default) as 1e-4 took 1-4 min a fit and barely moved. unused for lasso cells
floor_min_ratio <- function(alpha) ifelse(alpha == 0, 1e-4, 1e-2)

#' penalty and path floor of one cell. no alpha column means lasso
cell_floor <- function(cell) {
  alpha <- if (is.null(cell[["alpha"]])) 1 else cell[["alpha"]]
  list(alpha = alpha, max_imbalance = floor_max_imbalance,
       min_ratio = if (is.null(cell[["min_ratio"]])) floor_min_ratio(alpha)
                   else cell[["min_ratio"]])
}
