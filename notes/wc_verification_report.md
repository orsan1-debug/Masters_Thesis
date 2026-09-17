# Wong–Chan study: independent verification (text vs reference code vs our code vs sbw)

Date: 2026-09-16. Everything below was read from the sources listed in Part D; nothing is taken from earlier chat context.

Citation keys used throughout:

- **Sim.py** = `WongChanapproxBalSim.py` (68 lines); **Tune.py** = `WongChanapproxBalTune.py` (437 lines). Both files are byte-identical in the GitHub repo (yixinwang/minimal-weights-public, commit 3f586eb, 2026-08-17) and in `replication_code2.zip` from the Biometrika supplement (md5 a4bf86… / b1f5d7…). Line numbers apply to both copies.
- **arXiv** = arXiv:1705.00998v3 [stat.ME] 25 Apr 2019, PDF page numbers.
- **Supp** = `mabw_supp.pdf` (Biometrika Supplementary Material, 18 pp.), PDF page numbers, margin line numbers where printed.
- **Pub** = `WangZubizarreta2020.pdf` (Biometrika 107(1):93–105), printed page numbers.
- **dgp.R**, **est.R** (= `R/estimators_cv.R`), **sum.R** (= `analysis/cv/wz_replication/summarise_wc_sbw.R`), **plot.R** (= `analysis/cv/wz_replication/plot_wc_att.R`), **grids.R** (= `analysis/cv/wz_replication/wc_grids.R`), **qmd** (= `analysis/cv/wz_replication/wc_findings.qmd`), **run_path.R** etc. (= `runs/cv/wz_replication/run_wc_*.R`), all under `Thesis 2026/Masters_Thesis/`.
- **sbw** = CRAN package sbw version 1.2 (2025-10-15), files under `R/`. Python indices are 0-based; every `Z[:,j]` is translated to the 1-based label Z(j+1).

---

## PART A. Transcriptions

### A1. Appendix D.4, arXiv v3 (pp. 39–40; Figure 4 p. 41). Supp D.4 (pp. 15–17, margin ll. 897–921; Fig. 4 p. 18) has the same equations and sentences (checked sentence-by-sentence and visually).

Text as printed (arXiv p. 39):

> "It starts with a ten-dimensional multivariate standard Gaussian random vector Z = (Z1, …, Z10)ᵀ for each observation. Then it generates ten observed covariates X = (X1, …, X10)ᵀ, where
> X1 = exp(Z1/2), X2 = Z2/{1 + exp(Z1)}, X3 = (Z1Z3/25 + 0.6)³, X4 = (Z2 + Z4 + 20)², Xj = Zj, j = 5, …, 10.
> The propensity score model is pr(T = 1 | Z) = exp(−Z1 − 0.1Z4)/{1 + exp(−Z1 − 0.1Z4)}.
> The study considers two outcome regression models. Model A is Y = 210 + (1.5T − 0.5)(27.4Z1 + 13.7Z2 + 13.7Z3 + 13.7Z4) + ε, and model B is Y = Z1 Z2³ Z3² Z4 + Z4|Z1|^0.5 + ε, where ε ∼ N(0, 1)."

(arXiv p. 40) "We generate a dataset of size N = 5000 and study both the average treatment effect and the average treatment effect on the treated estimates. (We take the size of the bootstrap samples as 1/10 of the original sample size. We default to 10 bootstrap samples for covariate balance evaluation. We balance the first and second moments of the covariates.)"

| Item | Text |
|---|---|
| Latent Z | Z = (Z1,…,Z10)ᵀ, multivariate standard Gaussian |
| Observed X | X1 = exp(Z1/2); X2 = Z2/{1+exp(Z1)}; X3 = (Z1Z3/25+0.6)³; X4 = (Z2+Z4+20)²; Xj = Zj, j = 5,…,10 |
| Propensity | pr(T=1|Z) = expit(−Z1 − 0.1 Z4) |
| Outcome A | Y = 210 + (1.5T − 0.5)(27.4Z1 + 13.7Z2 + 13.7Z3 + 13.7Z4) + ε |
| Outcome B | Y = Z1 Z2³ Z3² Z4 + Z4 |Z1|^0.5 + ε |
| Noise | ε ∼ N(0,1) |
| n | N = 5000 |
| Number of datasets | **not stated** ("We generate a dataset") |
| Basis balanced | "first and second moments of the covariates" |
| Bootstrap | sample size = 1/10 of the original sample size; 10 bootstrap samples |
| Estimands | ATE and ATT |
| How RMSE truth is defined | **not stated** |
| Estimator form (HT vs Hájek) | **not stated** |
| Dispersion measures | Table 4 rows: Absolute Deviation, Variance, Negative Entropy |
| Any filtering of errors | **not stated** |
| Which balance metric selects δ | **not stated** in D.4; Algorithm 1 (Section 4.2) defines the statistic — see below |
| Table 4 caption (arXiv p. 40 / Supp p. 17) | "Table 4: Root mean squared error in the Wong-Chan study. Approximate balancing often produce similar-or-better quality estimates than exact balancing." |
| Figure 4 caption (arXiv p. 41 / Supp p. 18) | "Figure 4: Bootstrapped covariate balance C_S and mean squared error for different values of δ for the average treatment effect on the treated in the Wong and Chan study. Using C_S to select δ as in Algorithm 1 coincides with or neighbors the optimal δ with the smallest error, especially in estimating the average treatment effect. (The horizontal axis start from δ = 0. The vertical dotted line indicates δ = K^{−1/2}, where K is the number of covariates being balanced. We recommend not choosing δ's bigger than K^{−1/2} because they likely break the assumptions required by the asymptotics. )" Panels: (a) Average treatment effect on the treated, (b) Average treatment effect. Axis: x = δ (log, ticks 10⁻³…10⁰), y = "MSE / Cov. Bal." Legend: Abs. Dev. / Variance / Neg. Ent. × Cov. Bal. / MSE. |

**Algorithm 1, arXiv v3 p. 11** (verified visually):

> For each δ in a grid 𝒟 ⊂ [0, K^{−1/2}] of candidate imbalances
>   Compute {w_i}ⁿ_{i=1} by solving Problem (1)
>   For each k ∈ {1, …, K}
>     Draw a bootstrap sample 𝒦_k from the original data
>     Evaluate covariate balance C_k on the sample 𝒦_k,
>       C_k := ||{Σ_{i∈𝒦_k} w_i Z_i B_k(X_i)}/(Σ_{i∈𝒦_k} w_i Z_i) − Σ_{i=1}^n B_k(X_i)/n||₂ / sd{B_k(X)}
>   Compute the mean covariate balance, C_S(δ) := Σ_{k=1}^K C_k / K
> Output δ* = argmin_{δ∈𝒟} C_S(δ)

**Algorithm 1, Pub p. 101** (verified visually):

> For each δ in a grid 𝒟 ⊂ [0, K^{−1/2}] of candidate imbalances
>   Compute {w_i}ⁿ_{i=1} by solving problem (1)
>   For each j ∈ {1, …, J}
>     Draw a bootstrap sample 𝒮_j from the original data
>     Evaluate covariate balance C_j on the sample 𝒮_j,
>       C_j := Σ_{k=1}^K ||{Σ_{i∈𝒮_j} w_i Z_i B_k(X_i)}/(Σ_{i∈𝒮_j} w_i Z_i) − Σ_{i=1}^n B_k(X_i)/n||₂ / sd{B_k(X)}
>   Compute the mean covariate balance, C̄(δ) := Σ_{j=1}^J C_j / J
> Output δ* = argmin_{δ∈𝒟} C̄(δ)

Differences arXiv → Pub: loop index k∈{1..K} → j∈{1..J}; sample 𝒦_k → 𝒮_j; no Σ_k in front of the norm → Σ_{k=1}^K in front; C_S(δ) = Σ_k C_k/K → C̄(δ) = Σ_j C_j/J. Neither version states which units sd{B_k(X)} is taken over, the bootstrap sample size, or the number of draws (J or K). Both: weighted sum is over w_i Z_i (the weighted/treated group) in the bootstrap sample; the target is the full-sample mean Σ_{i=1}^n B_k(X_i)/n. Pub p. 101 text after the box: "We recommend choosing values of δ smaller than K^{−1/2}".

### A2. Sim.py, as executed

| Item | Code (line) | 1-based Z label |
|---|---|---|
| Z | `Z = npr.normal(size=(n, 10))` (L8) | Z1..Z10 iid N(0,1) |
| X1 | `X[:,0] = np.exp(Z[:,0]/2.)` (L10) | exp(Z1/2) |
| X2 | `X[:,1] = Z[:,1]/(1+np.exp(Z[:,0]))` (L11) | Z2/(1+exp(Z1)) |
| X3 | `X[:,2] = (Z[:,0]*Z[:,2]/25.+0.6)**3` (L12) | (Z1 Z3/25 + 0.6)³ |
| X4 | `X[:,3] = (Z[:,1]+Z[:,3]+20)**2` (L13) | (Z2 + Z4 + 20)² |
| X5..X10 | `X[:,4:] = Z[:,4:]` (L14) | Z5..Z10 |
| Propensity | `p = np.exp(-Z[:,1]-0.1*Z[:,4]) / (1.+np.exp(-Z[:,1]-0.1*Z[:,4]))` (L18) | **expit(−Z2 − 0.1 Z5)** |
| T | `T = npr.binomial(1, p)` (L19) | Bernoulli(p) |
| Y, model A | `210 + (1.5*T-0.5) * (27.4*Z[:,1]+13.7*Z[:,2]+13.7*Z[:,3]+13.7*Z[:,4]) + npr.normal(size=n)` (L23–25) | **210 + (1.5T − 0.5)(27.4 Z2 + 13.7 Z3 + 13.7 Z4 + 13.7 Z5) + ε** |
| Y1, model A | same with T := 1, `npr.normal(size=n)` drawn again (L26–28) | own noise draw |
| Y0, model A | same with T := 0, `npr.normal(size=n)` drawn again (L29–31) | own noise draw |
| Y, Y1, Y0, model B | `Z[:,1]*(Z[:,2]**3)*(Z[:,3]**2)*Z[:,4] + Z[:,4]*(np.abs(Z[:,1]))**0.5 + npr.normal(size=n)`, three separate noise draws (L35–40) | **Z2 Z3³ Z4² Z5 + Z5 |Z2|^0.5 + ε**, no T |
| Noise | `npr.normal(size=n)` = N(0,1); Y, Y1, Y0 each get an independent draw (L25, L28, L31, L36, L38, L40) | so Y ≠ T·Y1 + (1−T)·Y0 |
| n | `WongChanSimA(n=5000)`, `WongChanSimB(n=5000)` (L60, L63); function default n=200 (L43, L49) unused | 5000 |
| Datasets | `N = 100` (L57); loop `for i in range(N)` (L59) writes `sim_datasets/<i>WongChanSimA.csv` and `…B.csv` (L58, L62, L65) | 100 per model, A and B from separate draws of Z (L60 vs L63) |
| Column order in CSV | `np.column_stack([Z, p, X, T, Y, Y1, Y0])` (L61, L64) | cols 0–9 Z, 10 p, 11–20 X, 21 T, 22 Y, 23 Y1, 24 Y0 |
| RNG seed | none set | not reproducible |

### A3. Tune.py

**(i) Optimisation problem.** ATT (`approxBalATT`, L12–53): variable `w` of length n_ctrl (L18). Objectives: `"l1"`: minimize `cp.norm(w − 1/n_ctrl, 1)` (L21); `"l2"`: minimize `cp.sum_squares(w − 1/n_ctrl)` (L23); `"entropy"`: minimize `−cp.sum(cp.entr(w))` = Σ w log w (L25). Constraints: `cp.sum(w) == 1` (L30); `0 <= w` (L32); for each of the n_cov columns i: `X[T==0][:,i]*w − np.mean(X[T==1][:,i]) <= tol * X[:,i].std()` and the reverse (L33–39). Scale: `X[:,i].std()` = numpy std over **all n units, ddof=0** (L36, L39); target = treated mean (L35, L37). ATE (`approxBalATE`, L56–136): two separate problems, control weights `w_c` (L66–87) and treated weights `w_t` (L102–123), same objectives (L67–72, L103–108), sum-to-one (L77, L113), non-negativity (L79, L115); each arm's weighted mean constrained to within `tol * X[:,i].std()` of `np.mean(X[:,i])` (mean over all units) (L81–86, L117–122), again pooled ddof=0 std. Solver: `prob.solve()` with cvxpy default (L46, L93, L129). Failure handling: on `SolverError` or `w.value is None` the weight vector is set to `−1 * ones` (L47–51, L94–98, L130–134).

**(ii) Tolerance grid.** `tol_vals = np.array([0, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1])` (L213–215): 11 values including 0 and values above K^{−1/2} = 20^{−1/2} = 0.2236.

**(iii) Bootstrap balance.** ATT (`covBalATTBootstrapEval`, L138–164): if any weight equals the −1 sentinel, return `1e+16` (L142, L164). `n_subset = int(prop * n_ctrl)` with `prop=prop` passed from the main loop where `prop = 0.1` (L144, L218, L249–251; the signature default `prop=0.5` L138 is not used); `smps=10` draws (L138 default, not overridden); `subsamples = [npr.choice(n_ctrl, n_subset, replace=True) …]` — drawn **from the control units only, with replacement** (L147). Statistic per draw: `cov_dif = (w_c[sub]·X[T==0][sub] / Σ w_c[sub] − mean(X[T==1], axis=0)) / sd`, with `sd = np.std(X, axis=0)` over **all units, ddof=0** (L145, L149–150) — weighted control-subsample mean vs full treated mean. Metrics: `"l1"` = Σ_k |cov_dif_k| (L153); `"l2"` = Σ_k cov_dif_k² (L155); `"linf"` = max_k |cov_dif_k| (L157). Returned value: `np.nanmean(dif)` over the 10 draws (L161). ATE (`covBalATEBootstrapEval`, L166–202): combine `w[T==0] = w_c`, `w[T==1] = w_t` (L173–175); `n_subset = int(prop * n_all)` drawn **from all n units with replacement**, 10 draws (L176–177); per draw `cov_dif = (weighted mean of X over control units in the subsample − weighted mean of X over treated units in the subsample) / sd`, `sd = np.std(X, axis=0)` (L178, L182–188) — **treated-vs-control within the subsample, no population target**. Same three metrics (L190–195); `np.nanmean` (L199); `1e+16` if either arm failed (L170, L202).

**(iv) Point estimates.** ATT: `np.mean(Y[T==1]) − w_c.dot(Y[T==0])` (L257). ATE: `w_t.dot(Y[T==1]) − w_c.dot(Y[T==0])` (L274). Weights sum to one by constraint, so no separate normalisation is applied.

**(v) Truth.** `true_att = np.mean(Y1[T==1] − Y0[T==1])` (L244); `true_ate = np.mean(Y1 − Y0)` (L245), per dataset, from the stored Y1, Y0 columns (cols 23, 24; L237–239). Because Y1 and Y0 carry independent N(0,1) noise (A2), for model A this is the sample ATT/ATE of 1.5·g plus the mean noise difference; for model B it is the mean noise difference alone (no treatment term in B).

**(vi) From squared errors to reported RMSE.** Per objective, per dataset i, `att_sqerr = [(est − true_att)² for each of the 11 tol values]` (L257–259) → array N×11 (L279). `att_sqerrs[att_sqerrs > 1e2] = np.nan` (L284) and the same for the three covbal arrays (L285–287; this also removes the 1e16 failure sentinel). Same for ATE (L289–297). Per-tol mean squared error `nanmean(att_sqerrs, axis=0)` stored as `<obj>_err` (L311) and `nanmean(covbal)` as `<obj>_covbal_<metric>_mean` (L313–318). "Exact" RMSE = `sqrt(nanmean(att_sqerrs[:,0]))` — the tol = 0 column (L339). "Approx." RMSE, one per metric: for each dataset i, `loc = nanargmin(att_covbal_<metric>s[i,:])` over all 11 tol values (including 0 and the values above K^{−1/2}); take `att_sqerrs[i, loc]`; `sqrt(nanmean(...))` (L340–345). Three approx values are printed side by side with the exact value (L347). ATE identical (L349–357). Note `nanargmin` raises on an all-NaN row; no handling in the code. The Sim/Tune scripts share no seed, so the datasets in `sim_datasets/` are not reproducible from the code.

**(vii) Which metric feeds what.** Printed RMSE (L347, L357): exact, approx-l1, approx-l2, approx-linf — the code does not mark which of the three is the "Approx." column of Table 4. Figures (L365–437): plot `<obj>_covbal_l2_mean` (the **l2** metric: Σ_k standardised-difference²) as "Cov. Bal." and `<obj>_err` (MSE, not RMSE) as "MSE", against `tol` on log–log axes, for ATE (fig1, L400) and ATT (fig2, L437); values > 1e2 are set to inf before plotting (L362–363).

### A4. `dgp_wc()` in dgp.R (L221–248)

| Item | Code (line) |
|---|---|
| Z | `Z <- matrix(rnorm(n * 10), n, 10)` (L224): Z1..Z10 iid N(0,1) |
| Propensity | `eta <- -Z[, 1] - 0.1 * Z[, 4]` (L225); `e <- 1 / (1 + exp(-eta))` (L226) = expit(−Z1 − 0.1 Z4) |
| W | `rbinom(n, 1, e)` (L227) |
| g | `27.4 * Z[, 1] + 13.7 * Z[, 2] + 13.7 * Z[, 3] + 13.7 * Z[, 4]` (L229) |
| Outcome A | `210 + (1.5 * W - 0.5) * g` (L231) |
| Outcome B | `Z[, 1] * Z[, 2]^3 * Z[, 3]^2 * Z[, 4] + Z[, 4] * abs(Z[, 1])^0.5` (L232) |
| Noise | `eps <- rnorm(n)` once (L234); `Y <- vapply(outcome, \(o) f(o) + eps, numeric(n))` (L235): one N(0,1) draw shared by columns A and B; no Y1/Y0 columns |
| X | `X <- Z` (L237); if `misspec` (default TRUE, L221): `X[,1] <- exp(Z[,1]/2)`, `X[,2] <- Z[,2]/(1+exp(Z[,1]))`, `X[,3] <- (Z[,1]*Z[,3]/25 + 0.6)^3`, `X[,4] <- (Z[,2]+Z[,4]+20)^2` (L239–242); columns 5–10 stay Z5..Z10 |
| Individual effects | `tau_i <- 1.5 * g` (L245); returns `tau = 0`, `tau_i`, `e`, `eta`, `e0 = 1 - e` (L246–247) |
| n, reps | arguments of the caller (see A5) |

`dgp_wc_overlap()` (L268–296) is the same with `eta <- overlap * (-Z[,1] - 0.1*Z[,4])` (L273) and `eps <- sigma * rnorm(n)` (L282); defaults `overlap = 1, sigma = 1` (L268–269). `archive/cv_extension/dgp_wc_overlap.R` L1–20 is a copy of the same header and body.

### A5. Our sbw estimator and run

The sbw arm is defined in est.R; the script that produced `results/cv/wz_replication/wc_sbw_v2.csv.gz` is **not in the repository** (`notes/repo_state.md` L62: "wz producers missing: wc_att_v2, wc_sbw_v1/v2, wc_tuner_*"; `archive/cv_extension/run_wc.R` is 9 bytes, blank). Run-level facts below come from the CSV itself and from est.R defaults.

| Item | Code (line) |
|---|---|
| sbw call | `sbw(dat, ind = "W", bal = list(bal_cov = bal_cov, bal_alg = FALSE, bal_tol = delta, bal_std = "target"), sol = list(sol_nam = solver), par = list(par_est = estimand), mes = FALSE)` (est.R L22–26); `solver = "quadprog"` (L21, L52); `wei` not passed → sbw's formal default `wei = list(wei_sum = TRUE, wei_pos = TRUE)` (sbw.R L218); `out` not passed; `bal_gri`, `bal_sam` irrelevant with `bal_alg = FALSE` |
| Weights returned | `fit$dat_weights$sbw_weights` (L27), in data order |
| Basis | `X <- cbind(d$X, d$X^2)` (L54), 20 columns named x1..x10, x1sq..x10sq (L55); `bal_cov = colnames(X)` (L74) — first and second moments of the 10 observed covariates, K = 20 |
| Tolerance scale | `bal_std = "target"` (L24) → sbw multiplies `bal_tol` by `sd_target` (propar.R L52–53): ATT → sd (R `sd`, ddof=1) of the **treated** units (sbwcaufix.R L59); ATE → sd of **all** units, applied to each arm (sbwcaufix.R L42–43). `delta` is scalar and is replicated across the 20 columns (propar.R L45–48) |
| Grid | `grid = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)` (L50–51), 8 values; comment L49: "stays below K^(-1/2) = 0.22" |
| Exact arm | `delta_exact = 1e-4` (L52), fitted with `fit(delta_exact)` (L80); comment L11–13: quadprog rejects δ = 0 |
| Selection statistic `cstat()` | L33–45: for each arm `a` in `arms` (ATT: `0`; ATE: `c(0, 1)`, L71): `idx <- which(W == a)`; `m <- round(frac * length(idx))` with `frac = 0.1` (L33, L37) — **0.1 × arm size, from that arm only, with replacement** (L39); `B = 10` draws (L33, L38); `wm <- colSums(w[i] * B_mat[i, ]) / sum(w[i])` (L40); `total <- total + mean(abs(wm - target) / s)` (L41); return `total / B` (L44). `s <- apply(X, 2, sd)` = sd over **all units, ddof=1** (L57). Target: ATE → `colMeans(X)`; ATT → `colMeans(X[W == 1, ])` (L72). Statistic = mean over the 20 columns of the absolute standardised difference, summed over arms for ATE |
| δ selection | `cs <- vapply(w_grid, \(w) if (is.character(w)) Inf else cstat(...))` (L77–78); `k <- which.min(cs)` (L79); approx row uses `w_grid[[k]]`, `grid[k]` (L81) |
| Estimate | `tau_hat = colSums(w * W * d$Y) - colSums(w * (1 - W) * d$Y)` (L66) after `stopifnot(abs(sum(w * W) - 1) < 1e-6, abs(sum(w * (1 - W)) - 1) < 1e-6)` (L64); for ATT sbw sets every treated weight to 1/n₁ (sbwcaufix.R L66), so this is mean(Y_treated) − Σ w_c Y_c |
| Failure handling | `tryCatch(..., error = \(e) conditionMessage(e))` (L73–75) → row with `tau_hat = NA`, `error = message` (L60–63) |
| Output rows per rep | 2 estimands × 2 balance × outcome columns (L70–83) |
| Truth | not computed in est.R. CSV columns `sate`, `satt` (per rep) come from the missing run script. In sum.R: population truth `patt = -86.31 * c0` from Stein's lemma (L18–20) for ATT-A, 0 for ATE and for B (L24–25); sample truth `sate`/`satt` from the CSV, 0 for B (L26–27). In qmd L49–50 and L94, plot.R L26, L30: 0 for B, 0 for ATE, `satt` for ATT-A. `satt = mean(d$tau_i[W == 1])` in every run script (run_path.R L87; run_overlap.R L79; run_noise.R L79; run_overlap_n.R L81; run_basis.R L91) — noise-free |
| RMSE | sum.R L37 `sqrt(mean(err^2))`, MCSE L38 `sd(err^2)/sqrt(n())/(2*sqrt(mean(err^2)))`; qmd L52–54 same; no filtering, no `na.rm` |
| Reps, n | `wc_sbw_v2.csv.gz`: 1000 reps × 8 rows, 0 NA `tau_hat`, 0 error rows; approx `delta` ∈ {0.001, …, 0.2} (ATE never picked 0.2; ATT picked 0.2 in 44 of 2000 rows), exact `delta` = 1e-4 throughout — consistent with the est.R defaults. n is not stored in the CSV; qmd L27 states n = 5000; the fixed-n run scripts set `n <- 5000` (run_path.R L23, run_overlap.R L16, run_noise.R L15, run_basis.R L17; run_overlap_n.R varies n by cell) with `master_seed <- 20260903` (run_path.R L25), and run_path.R L15/L33 states the streams are identical to wc_sbw_v2. Not verifiable from the CSV alone. `wc_sbw_v1.csv.gz` has 10 reps |
| Comments citing paper settings | est.R L4–6 "B = 10 bootstrap samples of size n / 10, first and second moments, K = 20"; L7–9 "sbw::sbw(bal_alg = TRUE) bootstraps the full arm size bal_sam times and the subsample size cannot be changed, so the tuning loop is written out here around fixed-tolerance sbw fits" |

The balnet run scripts (run_path.R etc.) reuse `cstat()` with `arms = 0`, `target = colMeans(X[W == 1, ])`, `s = apply(X, 2, sd)` on ATT-scale balnet weights (run_path.R L75–79) and grids of 8 (run_path.R L26, run_overlap.R L20), 9 (run_overlap_n.R L21, run_basis.R L21: adds 0.5) or 10 values (run_noise.R L19: adds 0.5, 1). They do not call sbw.

### A6. sbw with `bal_alg = TRUE` (sbw 1.2)

| Item | Code (line) |
|---|---|
| Entry | `sbw()` sets `bal$bal_tol = bal$bal_gri` and calls `.sbwtun` (sbw.R L259–263) → `.sbwcautun` (sbwtun.R L13–14) |
| Grid default | `bal_gri = c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)` (sbw.R L218, L232–233); must be all positive (sbwauxtun.R L7) |
| Draws default | `bal_sam = 1000` (sbw.R L218, L238–239) |
| Target | ATE: `bal$bal_tar = colMeans(dat[, bal_cov])` over all units (sbwcautun.R L25–26); ATT: `colMeans` of treated (L27–28) |
| ATT flow | `sd_target = sd` of treated columns (L65); `.sbwauxtun(dat_level[[1]], …)` on the **control** data frame (L66); treated units get weight `1/(n − n₀)` (L72) |
| ATE flow | `sd_target = sd` over all units (L46); `.sbwauxtun` run on **each arm separately** with `run = FALSE` (L47); `cstat = colSums(rbind(cstat_control, cstat_treated))` — **summed over the two arms** per grid value (L48); `bal$bal_tol = bal$bal_gri[which.min(cstat)]` (L49); refit both arms at that tolerance (L50) |
| Per-arm tuning `.sbwauxtun` | for each `tol` in `bal_gri`: fit `.sbwauxfix` (L9–11); `n = length(weights)` = size of the arm being weighted (L13); `for (k in 1:bal$bal_sam)` (L15); `ind = sample(n, n, replace = TRUE)` — **resample of the full arm size, from that arm, with replacement** (L16); statistic per draw = `mean(abs((tᵀ(weights[ind]) %*% dat[ind, bal_cov] / sum(weights[ind]) − bal_tar) / <sd>))` — mean over covariates of the absolute standardised difference between the weighted resample mean and the target (L18, L20, L22); `<sd>` = `apply(dat[, bal_cov], 2, sd)` of the **weighted arm** when `bal_std = "group"` (L17–18), `sd_target` when `"target"` (L19–20), none when `"manual"` (L21–22); C_k accumulated and divided by `bal_sam` (L25); `bal_tol = bal_gri[which.min(Cstat.object)]` (L27–28) |
| Objective | `nor = "l_2"` (sbwauxfix.R L15); quadprog `Dmat = 2*diag(n)`, `c = 0` (sbwpri.R L217–220): minimise Σ wᵢ² subject to Σ wᵢ = 1 when `wei_sum = TRUE` and wᵢ ≥ 0 when `wei_pos = TRUE` (sbwauxfix.R L3–13; propar.R L108–112, L163–165) |
| Constraint scaling | `bal_tol` scalar → replicated (propar.R L45–48); `"group"`: `sd(bal_cov columns of the weighted arm) * bal_tol` (L50–51); `"target"`: `sd_target * bal_tol` (L52–53); bounds `bal_tar ± bal_tol` (L147–152) |
| δ = 0 | quadprog path stops if all tolerances are 0 (sbwauxfix.R L73–76) |

---

## PART B. Difference tables

Cells: **text** = A1 (arXiv D.4 / Algorithm 1, arXiv p. 11 unless noted); **code** = A2/A3; **ours** = A4/A5; **sbw** = A6 (`bal_alg = TRUE`, defaults). "n/a" = the source has no such item.

### B1. Data-generating process

| item | text | code | ours | sbw package | identical? |
|---|---|---|---|---|---|
| Z | 10-dim standard Gaussian (arXiv p. 39) | `npr.normal(size=(n,10))` (Sim L8) | `matrix(rnorm(n*10), n, 10)` (dgp.R L224) | n/a | text = code = ours |
| X1 | exp(Z1/2) (p. 39) | exp(Z1/2) (Sim L10) | exp(Z1/2) (L239) | n/a | yes |
| X2 | Z2/{1+exp(Z1)} | Z2/(1+exp(Z1)) (L11) | Z2/(1+exp(Z1)) (L240) | n/a | yes |
| X3 | (Z1Z3/25+0.6)³ | (Z1Z3/25+0.6)³ (L12) | (Z1Z3/25+0.6)³ (L241) | n/a | yes |
| X4 | (Z2+Z4+20)² | (Z2+Z4+20)² (L13) | (Z2+Z4+20)² (L242) | n/a | yes |
| X5..X10 | Zj, j=5..10 | Z5..Z10 (L14) | Z5..Z10 (L237, untouched) | n/a | yes |
| Propensity Z indices | −Z1 − 0.1 Z4 (p. 39) | **−Z2 − 0.1 Z5** (Sim L18) | −Z1 − 0.1 Z4 (dgp.R L225) | n/a | text = ours; **code differs** |
| Propensity link | expit | expit (L18) | expit (L226) | n/a | yes |
| Outcome A Z indices | 27.4Z1+13.7Z2+13.7Z3+13.7Z4 (p. 39) | **27.4Z2+13.7Z3+13.7Z4+13.7Z5** (Sim L24) | 27.4Z1+13.7Z2+13.7Z3+13.7Z4 (L229) | n/a | text = ours; **code differs** |
| Outcome A form | 210 + (1.5T−0.5)(·) + ε | same (L23–25) | `210 + (1.5*W − 0.5)*g` (L231) | n/a | yes |
| Outcome B Z indices | Z1 Z2³ Z3² Z4 + Z4|Z1|^0.5 | **Z2 Z3³ Z4² Z5 + Z5|Z2|^0.5** (Sim L35) | Z1 Z2³ Z3² Z4 + Z4|Z1|^0.5 (L232) | n/a | text = ours; **code differs** |
| Noise | ε ∼ N(0,1) | N(0,1); independent draws for Y, Y1, Y0 (Sim L25, L28, L31, L36–40) | N(0,1); one draw shared by A and B; no Y1/Y0 (L234–235) | n/a | distribution same; draw structure differs |
| n | N = 5000 (p. 40) | 5000 (Sim L60, L63) | 5000 (run_path.R L23, run_overlap.R L16, run_noise.R L15, run_basis.R L17; qmd L27); not stored in wc_sbw_v2 | n/a | yes as stated |
| reps | not stated | N = 100 (Sim L57; Tune L210) | 1000 (wc_sbw_v2.csv.gz; run scripts `N_SIM` default 1000) | n/a | **code 100 vs ours 1000**; text silent |
| A and B on same draw? | not stated | separate datasets (Sim L60 vs L63) | same Z, W, eps (L235) | n/a | differs |
| Seed | not stated | none | `master_seed = 20260903`, L'Ecuyer streams (run_path.R L25, L34–38); wc_sbw_v2 producer missing | n/a | — |
| Sample truth stored | not stated | Y1, Y0 columns (Sim L61) | `tau_i = 1.5*g` (L245); `satt = mean(tau_i[W==1])` in run scripts | n/a | differs (see B3 truth) |

### B2. Balancing problem and tuning

| item | text | code | ours | sbw package | identical? |
|---|---|---|---|---|---|
| Basis | first and second moments (p. 40) | `[cov, cov**2]`, 20 cols (Tune L242) | `cbind(d$X, d$X^2)`, 20 cols (est.R L54) | user-supplied `bal_cov` | yes (K = 20) |
| Objective ("Variance" row) | f(w) = (w − 1/r)² (Pub p. 101–102, Sec. 4.3) | `sum_squares(w − 1/n_ctrl)` (Tune L23) | sbw objective (below) | Σ wᵢ² (sbwpri.R L220) | code and sbw differ by the constant −1/n given Σw = 1: Σ(w−1/n)² = Σw² − 1/n; same minimiser under Σw = 1 |
| Other objectives | abs. dev., neg. entropy | `norm(w−1/n,1)` (L21), `−Σ entr(w)` (L25) | not run | quadprog: l_2 only (sbwpri.R L197–199) | ours/sbw: Variance only |
| Weight constraints | Σw = 1, w ≥ 0 optional (Pub p. 96) | Σw = 1 (L30), w ≥ 0 (L32) | sbw defaults `wei_sum = TRUE, wei_pos = TRUE` (sbw.R L218) | same | yes |
| Tolerance scale, ATT | Algorithm 1: sd{B_k(X)}, units not stated | `tol * X[:,i].std()`: **all units, ddof=0** (Tune L36, L39) | `bal_std = "target"` → **treated** units, ddof=1 (est.R L24; sbwcaufix.R L59; propar.R L52–53) | default `"group"` → **control** units, ddof=1 (propar.R L50–51) | all three differ |
| Tolerance scale, ATE | as above | `tol * X[:,i].std()`: all units, ddof=0 (L83, L86, L119, L122) | `"target"` → all units, ddof=1 (sbwcaufix.R L42) | `"group"` → each arm's own sd | code ≈ ours up to ddof; sbw default differs |
| Balance target, ATT | Σᵢ B_k(Xᵢ)/n (full sample) in Alg. 1; D.4 silent | treated mean (Tune L35) | treated mean (sbw `par_est = "att"`; sbwcaufix.R L58) | treated mean (sbwcautun.R L28) | code = ours = sbw; Alg. 1 as printed uses the full-sample mean |
| Balance target, ATE | full-sample mean | full-sample mean, each arm (L82, L118) | full-sample mean, each arm (sbwcaufix.R L41) | full-sample mean (sbwcautun.R L26) | yes |
| Grid | 𝒟 ⊂ [0, K^{−1/2}] (p. 11); recommend δ < K^{−1/2} | {0, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1} (Tune L213–215): includes 0 and 3 values > 0.2236 | approx: {1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2} (est.R L50–51); exact: 1e-4 (L52) | {1e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1} (sbw.R L233) | all differ |
| Exact arm | δ = 0 (Table 4 "Exact") | tol = 0 column (Tune L339, L349) | δ = 1e-4 (est.R L52) | 1e-4 is the smallest grid value; 0 rejected by quadprog (sbwauxfix.R L74–75) | differs |
| Resample source, ATT | "from the original data" (p. 11) | controls only (Tune L147) | controls only (est.R L36, L71) | controls only (sbwcautun.R L66; sbwauxtun.R L16) | code = ours = sbw |
| Resample source, ATE | "from the original data" | **all units jointly** (Tune L177) | each arm separately (est.R L35–36) | each arm separately (sbwcautun.R L47; sbwauxtun.R L16) | code differs from ours/sbw |
| Resample size | 1/10 of the original sample size (p. 40) | ATT `int(0.1·n_ctrl)` (L144); ATE `int(0.1·n_all)` (L176) | `round(0.1·arm size)` per arm (est.R L37) | full arm size n (sbwauxtun.R L16) | text says 1/10 of original sample; code: 1/10 of controls (ATT) / 1/10 of all (ATE); ours: 1/10 of each arm; sbw: full arm |
| With/without replacement | not stated | with (L147, L177) | with (est.R L39) | with (L16) | code = ours = sbw |
| Draws | 10 (p. 40); Alg. 1 K (arXiv) / J (Pub), unspecified | 10 (Tune L138, L166 defaults) | 10 (est.R L33) | 1000 (sbw.R L239) | text = code = ours; sbw differs |
| Bootstrap sd | sd{B_k(X)} | `np.std(X, axis=0)`, all units, ddof=0 (L145, L178) | `apply(X, 2, sd)`, all units, ddof=1 (est.R L57) | `"group"`: weighted arm; `"target"`: sd_target (sbwauxtun.R L17–20) | code ≈ ours (ddof); sbw differs |
| Balance statistic | arXiv: ||·||₂/sd per k, averaged over draws; Pub: Σ_k ||·||₂/sd{B_k}, averaged over draws | l1 = Σ_k|d_k|, l2 = Σ_k d_k², linf = max_k|d_k| (L152–157, L190–195); figures use l2; table metric not marked | mean_k |d_k| (est.R L41) | mean_k |d_k| (sbwauxtun.R L18, L20) | ours = sbw (ℓ1/K); code has three, none is ℓ1/K; text ℓ2 |
| ATE criterion | Alg. 1: one weighted group vs full-sample mean (w_i Z_i) | weighted **control vs treated** subsample means (L186–187) | each arm vs full-sample mean, **summed** over arms (est.R L34–44) | each arm vs full-sample mean, `colSums` over arms (sbwcautun.R L48) | ours = sbw; code differs |
| Failure value | n/a | 1e16 (L164, L202), later NaN (L285–287) | Inf (est.R L77) | n/a (errors propagate) | — |
| Selection rule | argmin over 𝒟 of mean balance | per dataset `nanargmin` over all 11 tol (L341, L343, L345) | `which.min(cs)` over 8 (est.R L79) | `which.min` over grid (sbwauxtun.R L28; sbwcautun.R L49) | same rule, different grids/statistics |
| Solver | not stated | cvxpy default (L46) | quadprog (est.R L21, L52) | quadprog default (sbw.R L218) | — |

### B3. Estimator, truth, error processing

| item | text | code | ours | sbw package | identical? |
|---|---|---|---|---|---|
| ATT estimator | not stated | mean(Y_t) − w_cᵀY_c (Tune L257) | Σ w W Y − Σ w(1−W)Y with treated w = 1/n₁ (est.R L66; sbwcaufix.R L66) = mean(Y_t) − w_cᵀY_c | weights only; `estimate()` not read | code = ours |
| ATE estimator | not stated | w_tᵀY_t − w_cᵀY_c (L274) | same form (est.R L66) | — | yes |
| Normalisation | not stated | none beyond Σw = 1 constraint | `stopifnot` both arms sum to 1 (est.R L64) | Σw = 1 constraint | HT = Hájek in both |
| Truth, ATT | not stated | `mean(Y1[T==1] − Y0[T==1])`, includes noise difference (L244) | `satt = mean(tau_i[W==1])`, noise-free (run scripts); population `patt` in sum.R L20 | n/a | differs |
| Truth, ATE | not stated | `mean(Y1 − Y0)` (L245) | 0 (qmd L49; sum.R L24 population) or `sate` (sum.R L26; CSV column) | n/a | differs |
| Truth, model B | not stated | mean of noise difference (Y1−Y0 has no treatment term) | 0 (qmd L49; sum.R L24, L26) | n/a | differs |
| Error filtering | not stated | squared errors > 100 → NaN (L284, L294); covbal > 100 → NaN (L285–287, L295–297); nanmean | none (sum.R L37; qmd L52); 0 NA in wc_sbw_v2 | n/a | differs |
| RMSE | "Root mean squared error" (Table 4 caption) | `sqrt(nanmean(sqerr))` (L339–345) | `sqrt(mean(err^2))` (sum.R L37; qmd L52) | n/a | same formula modulo NaN |
| Datasets averaged | not stated | ≤ 100 | 1000 | n/a | differs |
| MCSE | not reported | none | `sd(err²)/√n/(2·RMSE)` (sum.R L38; qmd L54; grids.R L27) | n/a | — |

---

## PART C. Table 4 check

**Table 4 (arXiv p. 40 = Supp p. 17), copied; bold as printed:**

(a) Average treatment effect on the treated

| Minimize | A Exact | A Approx. | B Exact | B Approx. |
|---|---|---|---|---|
| Absolute Deviation | 0.67 | **0.66** | **0.26** | **0.26** |
| Variance | **0.72** | 0.79 | 0.26 | **0.25** |
| Negative Entropy | **0.78** | 0.89 | **0.25** | **0.25** |

(b) Average treatment effect

| Minimize | A Exact | A Approx. | B Exact | B Approx. |
|---|---|---|---|---|
| Absolute Deviation | 0.47 | **0.45** | **0.23** | 0.24 |
| Variance | 1.35 | **0.51** | 0.31 | **0.21** |
| Negative Entropy | **0.44** | 0.52 | **0.21** | **0.21** |

Caption: "Table 4: Root mean squared error in the Wong-Chan study. Approximate balancing often produce similar-or-better quality estimates than exact balancing."

**Hard-coded paper values, qmd L42–47 (`paper` tribble):**

| qmd row | qmd value | Table 4 cell | Table 4 value | match |
|---|---|---|---|---|
| ate, A, exact | 1.35 | (b) Variance, A Exact | 1.35 | yes |
| ate, A, approx | 0.51 | (b) Variance, A Approx. | 0.51 | yes |
| ate, B, exact | 0.31 | (b) Variance, B Exact | 0.31 | yes |
| ate, B, approx | 0.21 | (b) Variance, B Approx. | 0.21 | yes |
| att, A, exact | 0.72 | (a) Variance, A Exact | 0.72 | yes |
| att, A, approx | 0.79 | (a) Variance, A Approx. | 0.79 | yes |
| att, B, exact | 0.26 | (a) Variance, B Exact | 0.26 | yes |
| att, B, approx | 0.25 | (a) Variance, B Approx. | 0.25 | yes |

**Captions and text citing paper values:**

| location | statement | Table 4 | flag |
|---|---|---|---|
| qmd L143 caption template + L146 `show_model("A", "0.25 exact, 0.25 approx")` | "Model A, ATT … Their Neg. Ent. row: 0.25 exact, 0.25 approx." | (a) Negative Entropy, model A: **0.78 exact, 0.89 approx** | **mismatch** — 0.25/0.25 is the (a) Negative Entropy model **B** row |
| qmd L149 `show_model("B", "0.25 exact, 0.25 approx")` | "Model B, ATT … Their Neg. Ent. row: 0.25 exact, 0.25 approx." | (a) Negative Entropy, model B: 0.25, 0.25 | matches |
| qmd L72 | "Paper: 10% worse" (model A, ATT) | (a) Variance A: 0.72 → 0.79 (+9.7%) | consistent with the Variance row |
| qmd L27 | "n = 5000, 10 covariates, first and second moments balanced" | D.4 p. 39–40 | matches; "1000 paired reps" is ours, not the paper's |
| qmd L68 | "model B and both exact arms match Table 4 within 0.35; model A approximate does not (2.4 to 2.7 times theirs)" | from wc_sbw_v2 (recomputed: ATT-A approx 2.11 vs 0.79 = 2.7×; ATE-A approx 1.23 vs 0.51 = 2.4×; ATT-A exact 0.92 vs 0.72; ATE-A exact 1.00 vs 1.35; B arms 0.23–0.29 vs 0.21–0.31) | arithmetic consistent with the CSV and Table 4 |
| qmd L170 | "sbw, their selector (λ = 0.04 on sbw's scale, RMSE = 0.23)" | ours; 0.04 equals the **mean** approx δ for ATT in the CSV (0.0397); the **median** is 0.02, which is what qmd L68 quotes | internal inconsistency between L68 (median 0.02) and L170 (0.04), not a paper-value issue |
| sum.R L2, est.R L2 | "Table 4 'Variance' rows" | rows exist as copied | consistent |
| plot.R, grids.R | no paper values hard-coded | — | — |

---

## PART D. What was and was not verified

Verified from the source itself:

- Sim.py and Tune.py: read in full from the GitHub clone and from the supplement zip; identical.
- arXiv v3: PDF obtained (601,677 bytes, pdfTeX, dated 2019-04-26) via your desktop's browser pane (the cloud proxy blocks arxiv.org); Algorithm 1 (p. 11), D.4 (pp. 39–40), Table 4 (p. 40) and Figure 4 (p. 41) read from `pdftotext` output and checked visually on page renders.
- Published article: `readings/WangZubizarreta2020.pdf` (the OUP page itself is paywalled); Section 4.2 / Algorithm 1 (p. 101) read from text and page render. The project doc `WangZubizarreta2020.md` was not used for any cell.
- Supplement: `readings/asz050-suppl_data/mabw_supp.pdf` (Quartz, 2019-09-20); D.4 pp. 15–17, Fig. 4 p. 18, read from text and page render.
- sbw: CRAN mirror on GitHub, version 1.2; sbw.R, sbwtun.R, sbwcautun.R, sbwauxtun.R, sbwauxfix.R read in full, plus sbwcaufix.R, propar.R and the quadprog block of sbwpri.R. The rdrr.io pages themselves were not fetched. The sbw version installed on your machine was not checked (`R/packages.R` L9 loads it without a version).
- Our files: dgp.R, estimators_cv.R, summarise_wc_sbw.R, plot_wc_att.R, wc_grids.R, wc_findings.qmd, the five run_wc_*.R, archive/cv_extension/dgp_wc_overlap.R and run_wc.R, notes/repo_state.md, wc_sbw_v1/v2 CSVs — all read in full.

Not verifiable:

- The script that produced `wc_sbw_v2.csv.gz` (and v1) is absent from the repository; n, seed and the exact call to `estimate_wc_sbw()` for that batch are inferred only from the CSV contents and est.R defaults, as marked in A5.
- `wc_sbw_v1_meta.rds` / `wc_sbw_v2_meta.rds` could not be opened here (no R runtime in the sandbox; the Python reader install timed out).
- Tune.py does not identify which of its three "Approx." values (l1, l2, linf) is the one in Table 4, and D.4 does not say.
- The `sim_datasets/` CSVs behind Table 4 are not in the repo; Table 4 cannot be regenerated from the code without them (no seed).
- sbw `estimate()` (`R/estimate.R`) was not read; our code does not call it.
