# *** Setup ***
# Validity tests for dgp2(); sources the live generator, never a copy.
root <- getwd()   # walk up to the project folder; runs from root or any subfolder on any machine
while (basename(root) != "Masters_Thesis" && dirname(root) != root) root <- dirname(root)
stopifnot(basename(root) == "Masters_Thesis")
setwd(root)
source("R/dgp.R")

ess      <- function(w) sum(w)^2 / sum(w^2)
true_w   <- function(d) list(w1 = 1 / d$e[d$W == 1], w0 = 1 / d$e0[d$W == 0])   # e0 stays finite near e = 1
hajek    <- function(w, e, e0, y) {
  sum(w * y / e) / sum(w / e) - sum((1 - w) * y / e0) / sum((1 - w) / e0)
}
hajek_se <- function(w, e, e0, y) {   # e known; Kish var(Y) / ESS understates it (big weights sit where |y - mu| is big)
  arm <- function(wa, ya) sum(wa^2 * (ya - sum(wa * ya) / sum(wa))^2) / sum(wa)^2
  sqrt(arm(1 / e[w == 1], y[w == 1]) + arm(1 / e0[w == 0], y[w == 0]))
}

### ---- strength AND mu0 TOUCH ONLY Y ---- ###
# Hajek ignores mu0, HT does not (Ben-Michael et al. 2021, p. 22); 210 = Kang & Schafer 2007 level.
set.seed(3); y0 <- dgp2(1e4, p = 4, strength = 0)$Y
set.seed(3); d1 <- dgp2(1e4, p = 4)
set.seed(3); y2 <- dgp2(1e4, p = 4, strength = 2)$Y
set.seed(3); d3 <- dgp2(1e4, p = 4, mu0 = 210)
all.equal(y2 - d1$Y, d1$Y - y0)                                    # TRUE
range(d3$Y - d1$Y)                                                 # 210 210
identical(d3[c("W", "X", "e")], d1[c("W", "X", "e")])              # TRUE
hajek(d3$W, d3$e, d3$e0, d3$Y) - hajek(d1$W, d1$e, d1$e0, d1$Y)    # ~0

### ---- OVERLAP LEVELS: ONE DRAW PER LEVEL, n = 1e6 ---- ###
# p = 4: eta and Y use only the first s = 4 columns; p = 60 costs 480 MB per matrix.
levels <- c(0.25, 0.5, 1, 2)
dd <- lapply(levels, function(o) { set.seed(1); dgp2(1e6, p = 4, overlap = o, tau = 1) })

### ---- OVERLAP LEVELS: PROPENSITY RANGE ---- ###
# eta ~ N(0, 1 / overlap^2), so P(e outside [0.1, 0.9]) = 2 pnorm(-qlogis(0.9) overlap).
ps_range <- t(vapply(dd, function(d) c(
  e_q01   = unname(quantile(d$e, 0.01)),
  e_q99   = unname(quantile(d$e, 0.99)),
  outside = mean(d$e < 0.1 | d$e > 0.9)), numeric(3)))
signif(cbind(overlap = levels, ps_range, outside_th = 2 * pnorm(-qlogis(0.9) * levels)), 3)   # outside = outside_th

### ---- OVERLAP LEVELS: ESS SHARE AND MAX TRUE WEIGHT ---- ###
# Large-sample share per arm 1 / (1 + exp(1 / (2 overlap^2))) (lognormal mean, derived).
# At 0.25 single weights > 1e4 set the ESS (Austin & Stuart 2015, p. 3663): not a grid value.
ess_tab <- t(vapply(dd, function(d) {
  w <- true_w(d)
  c(ess1_share = ess(w$w1) / length(d$W), ess0_share = ess(w$w0) / length(d$W),
    w_max = max(w$w1, w$w0))
}, numeric(3)))
signif(cbind(overlap = levels, ess_tab, ess_th = 1 / (1 + exp(1 / (2 * levels^2)))), 3)   # = ess_th for overlap >= 0.5

### ---- ORACLE RECOVERS tau AT EVERY OVERLAP LEVEL ---- ###
# First pilot check (Morris et al. 2019, sec. 5); a fail here is the DGP, not an estimator.
oracle <- t(vapply(dd, function(d) c(
  err = hajek(d$W, d$e, d$e0, d$Y) - d$tau,
  se  = hajek_se(d$W, d$e, d$e0, d$Y)), numeric(2)))
cbind(overlap = levels, signif(oracle, 3), pass = abs(oracle[, "err"]) < 3 * oracle[, "se"])   # pass = 1; one draw, so 3 se not an MCSE


### ---- misspec TOUCHES ONLY X (fix a: PS-wrong / OR-correct cell) ---- ###
# Kang & Schafer 2007 design: the analyst sees transformed covariates; W, e, Y unchanged.
# Oracle uses true e, so it must still pass; a fail means misspec leaked into e or Y.
set.seed(1); d_cor <- dgp2(1e6, p = 4, tau = 1)
set.seed(1); d_mis <- dgp2(1e6, p = 4, tau = 1, misspec = TRUE)
identical(d_mis[c("W", "e", "e0", "Y")], d_cor[c("W", "e", "e0", "Y")])   # TRUE
any(d_mis$X != d_cor$X)                                                    # TRUE
c(err = hajek(d_mis$W, d_mis$e, d_mis$e0, d_mis$Y) - d_mis$tau,
  se  = hajek_se(d_mis$W, d_mis$e, d_mis$e0, d_mis$Y))                     # |err| < 3 se



### ---- mu0 MOVES HT, NOT HAJEK ---- ###
# HT shift = mu0 * (mean(W / e) - mean((1 - W) / e0)): weight-sum imbalance times the level.
# Imai & Ratkovic 2014 Table 1, p. 11 is the motivating case (HT RMSE 2371 vs Hajek 12.7).
ht <- function(w, e, e0, y) mean(w * y / e) - mean((1 - w) * y / e0)
c(ht_shift    = ht(d3$W, d3$e, d3$e0, d3$Y) - ht(d1$W, d1$e, d1$e0, d1$Y),
  ht_th       = 210 * (mean(d1$W / d1$e) - mean((1 - d1$W) / d1$e0)),
  hajek_shift = hajek(d3$W, d3$e, d3$e0, d3$Y) - hajek(d1$W, d1$e, d1$e0, d1$Y))   # ht_shift = ht_th; hajek ~ 0