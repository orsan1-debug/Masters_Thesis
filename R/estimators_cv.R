# Algorithm 1 criterion from Wang and Zubizarreta 2020 for the Wong and Chan runs
# defaults are their settings (B = 10 bootstrap subsamples of size n / 10)

#' mean over B subsamples (frac x arm size) of the mean std abs imbalance of the weighted arms
#' vs the target means. summed over arms as sbw does for the ATE
cstat <- function(w, B_mat, W, arms, target, s, B = 10, frac = 0.1) {
  total <- 0
  for (a in arms) {
    idx <- which(W == a)
    m   <- round(frac * length(idx))
    for (b in seq_len(B)) {
      i  <- sample(idx, m, replace = TRUE)
      wm <- colSums(w[i] * B_mat[i, , drop = FALSE]) / sum(w[i])
      total <- total + mean(abs(wm - target) / s)
    }
  }
  total / B
}
