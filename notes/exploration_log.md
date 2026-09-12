# Exploration log

State of the simulation study as of 2026-09-12 (repo at `868f77f`). Compiled from the code in
`R/`, the run scripts in `runs/`, the results on disk, the old run ledger (47 rows, last committed
in `c88e082:registry.csv`) and the prose in the analysis documents. Findings are quoted verbatim
from those documents with file and line; nothing has been re-interpreted.

## 1. DGPs (`R/dgp.R`)

| function | signature | used by |
|---|---|---|
| `dgp1` | `dgp1(n, p, outcome = "linear", misspec = FALSE, covcor = "iid", overlap = 1)` | ipw, legacy (E1 batch `correctspecDGP1`) |
| `dgp2` | `dgp2(n, p, s = 4, outcome = "linear", misspec = FALSE, covcor = c("iid", "ar1"), overlap = 1, signs = c("pos", "neg", "mixed", "ks"), decay_ps = 1, decay_out = 1, outcome_set = c("fixed4", "track_s"), idx_ps = NULL, idx_out = NULL, treat_prop = 0.5, tau = 0, strength = 1)` | ipw, current (all other ipw batches) |
| `dgp_wc` | `dgp_wc(n, outcome = c("A", "B"), misspec = TRUE)` | cv / wz_replication (`run_wc_path.R`, sbw replication) |
| `dgp_wc_overlap` | `dgp_wc_overlap(n, overlap = 1, sigma = 1, outcome = c("A", "B"), misspec = TRUE)` | cv / wz_replication (overlap, noise, n, basis grids) |
| `gen_data` | `gen_data(n = 1000, p = 100, rho = 0.5, s_prop = 5, s_y = 5, overlap = c("bad", "moderate", "good"), sigma_y = 1, seed = NULL)` | cv / basic_dgp (every family; scripts wrap it in a local `gen_cell()` for "awful" overlap, numeric `c_prop`, density and spread variants) |
| `run_par` | `run_par(X, FUN, ...)` | cv / basic_dgp parallel helper (needs a global `cl`) |

- `dgp1`: logistic propensity on the first four of p iid or AR(1) covariates; outcomes linear /
  quad1 / exp on the same four columns; Kang-Schafer transforms when `misspec = TRUE`; tau = 0.
- `dgp2`: sparse propensity on s covariates with j^(-decay_ps) coefficients under a sign pattern,
  index normalised so sd(eta) = 1 / overlap, intercept solved for `treat_prop`; outcome loads on
  the first four propensity covariates ("fixed4") or all s ("track_s"); `strength` scales the
  signal; returns e and e0 separately. Validity tests: `tests/dgp2_validity.R`.
- `dgp_wc` / `dgp_wc_overlap`: Wong and Chan (2018) design as run by Wang and Zubizarreta (2020):
  10 latent normals, outcome models A (confounded, heterogeneous effect) and B (no effect),
  transformed observed covariates; the overlap variant scales the treatment logit and the noise SD.
- `gen_data`: p AR(1) covariates, s_prop propensity and s_y outcome coordinates overlapping on the
  first min(s_prop, s_y); target mu1 = E[Y(1)] = 0; overlap via the propensity scale
  (good 0.7, moderate 1.5, bad 2.5; "awful" = 4 and numeric `c_prop` set in the scripts).

## 2. Estimators

### ipw (`R/estimators_ipw.R`, one `estimate_all()` call per replication)

| row | what |
|---|---|
| `balnetcv` | balnet weights at the `cv.balnet` selected lambda (5 folds, `max.imbalance = 1e-4`) |
| `balnet0`, `balnet05`, `balnet10` | balnet weights at fixed lambda 0 (path endpoint), 0.05, 0.10 |
| `balnetrate` | balnet weights at lambda = sqrt(log(p) / n) (Wager 2024, s.7.2) |
| `glmnetcv_ht`, `glmnetcv_hajek` | `cv.glmnet` logistic propensity, Horvitz-Thompson and Hajek |
| `oracle_ht`, `oracle_hajek` | true propensity, Horvitz-Thompson and Hajek |
| diagnostics | `lam_balcv1/0`, `lam_end1/0`, `lam_glmcv`, `trunc05`, `nnz_*`, `smd1/0_*`, `prev`, `emin`, `emax`, `nout05`, `nout01`, `cvloss_*` |

Helpers: `ate_ht()`, `ate_hajek()`, `ate_bal()`, `max_smd()`, `cv_loss()`.

### cv / basic_dgp (selectors on the balnet path, treated arm, target mu1)

`cv.bloss` (cv.balnet, balance loss), `cv.smd` (cv.balnet, mean |SMD|), `cv.inf` (cv.balnet,
max |SMD|), `boot.smd`, `boot.inf` (cv.boot.balnet, 500 half-samples); reference: the RMSE along
the full lambda path (`est_path`) and its floor.

### cv / wz_replication (`R/estimators_cv.R` plus the run scripts)

The five selectors above on the ATT (target = "control"), plus `alg1` (Wang and Zubizarreta
Algorithm 1 via `cstat()`), `glmnet` and `glm` propensity ATT (Hajek on controls), `naive`
difference in means, and the sbw replication `fit_sbw()` / `estimate_wc_sbw()` (exact and
Algorithm-1-tuned stable balancing weights, quadprog solver).

## 3. Batches

### ipw (`runs/ipw/`, output `results/ipw/<batch>.csv.gz`, 1000 reps per cell unless noted)

| batch file | script | design | seed | ledger (old registry) |
|---|---|---|---|---|
| correctspecDGP1 | correctspecDGP1.R | dgp1, p 50, n 500..50k, outcomes linear/quad1/exp | 101 | rows 1-2, 7 |
| correctspecDGP2 | correct_misspec_dgp2.R | dgp2 correct, p 50, n 500..50k | 202 | rows 3, 18-19, 31; file on disk is the 7-cell (n to 200k), 2-outcome batch from `correctspec_dgp2_n800k.R`, balnet 0.0.3 (see phase4 report) |
| misspecDGP2 | correct_misspec_dgp2.R | dgp2 misspec, p 50, n 500..50k | 202 | rows 4, 20 |
| correctspecDGP2_ks / _pos / _mixed | (script not in repo) | dgp2 sign patterns | 202 | rows 8-10 |
| correctspecDGP2cvfix, misspecDGP2cvfix | cvfix.R | seed-shared reruns under the runtime cv.balnet rowMeans patch | 202 | rows 15-16 |
| overlapDGP2, overlapmisspecDGP2 | overlap_misspec_dgp2.R | n 500..10k x overlap 0.25..2 | 203 | rows 5-6, 11-12, 21-22 |
| dimensionDGP2, dimensionmisspecDGP2 | DimensionDGP2.R | n 500, 1000 x p/n 0.1..2 | 204 | rows 13-14, 23-24 |
| PoverMDGP2 | PoverMDGP2.R | n 500, 1000 x p/n 4, 8 | 205 | rows 17, 25-26 |
| sparsityDGP2 | Sparsity.R | n x s 4..50 x decay x covcor | 206 | rows 27-28 |
| instrumentsDGP2 | weak_instruments.R | n x s_instr 0..16 x decay | 207 | rows 29-30, 33-34 |
| correctspecDGP400k, correctspecDGP800k | correctspec_dgp2_n800k.R | dgp2 correct, n 400k / 800k | 203 / 204 | rows 35-37 |
| snrDGP2correct | SNR.R | n 1k, 10k x strength 0.25..4 | 203 | row 38 |
| prevDGP2correct | treat_prop.R | n 1k, 10k x treat_prop 0.05..0.5 | 204 | row 39 |
| sparsityE6DGP2, sparsityE6ar1DGP2 | sparsity_2.R | n 1000 x p 50..500 x s 4..256 x decay, iid / ar1 | 208 / 209 | rows 45, 47 |

### cv / basic_dgp (`runs/cv/basic_dgp/`, per-cell `results/cv/basic_dgp/<family>_<ii>_r<reps>.rds`)

| family | cells on disk | reps | script | grid |
|---|---|---|---|---|
| tune | 12 | 500 | run_noise_x_overlap_density_diluted__tune.R | overlap x s_y x sigma_y, lasso (results are the 10 Sep reruns; 2 Sep originals archived) |
| tune2 | 12 | 500 | run_noise_x_overlap_density_fixed__tune2.R | overlap x density x sigma_y, fixed confounding |
| tune4 (+ tune4_19 at 1000) | 21 | 500 / 1000 | run_noise_x_overlap_lasso__tune4.R, run_bad_overlap_noise10_1000reps__tune4_19.R | overlap (good..awful) x sigma_y 1..10, lasso |
| snr_enet | 20 | 500 | run_noise_x_overlap_enet__snr_enet.R | same grid, elastic net alpha 0.5 |
| tunea | 18 | 50 | run_alpha_x_overlap_x_density_50reps__tunea.R | alpha x overlap x density |
| tunea3 | 16 | 500 | run_ridge_enet_x_overlap_x_density__tunea3.R | ridge / EN x overlap x density |
| tunea4 | 8 | 200 | run_ridge_enet_lowered_floor__tunea4.R | ridge / EN, lowered path floor |
| alpha_bad | 4 | 200 | run_alpha_sweep_bad_overlap__alpha_bad.R | alpha sweep, bad overlap, sigma_y on shared fits |
| n_enet | 6 | 200 | run_n_x_overlap_enet__n_enet.R | n x overlap, elastic net |
| tunep2 | 18 | 50 | run_p_x_overlap_x_density__tunep2.R | p 500..2000 x overlap x density (smoke count) |
| spread | 36 | 200 | run_confounder_spread_x_penalty__spread.R | confounder spread s x overlap x alpha |
| dima | 8 | 200 | run_p_x_alpha_x_overlap__dima.R | p 100, 500 x alpha x overlap |
| dimhi | 12 | 200 | run_high_p_x_alpha_x_overlap__dimhi.R | p 1000, 2000 x alpha x overlap |
| snr_mb | 4 | 500 | run_noise_axis_moderate_bad__snr_mb.R | noise axis, moderate / bad overlap |
| ov_s1, ov_s1_m5, ov_1k | 6, 2, 3 | 500 / 500 / 1000 | run_overlap_axis__ov_s1.R, run_overlap_axis_maxit1e5__ov_s1_m5.R, run_overlap_axis_1000reps__ov_1k.R | overlap axis at sigma_y 1 |
| tunen2 | 16 | 500 | no script in repo | n x overlap x SNR, lasso |
| dimsnr | 8 | 200 | no script in repo | SNR at high p x alpha x overlap |
| snr_good | 7 | 500 | no script in repo | noise axis, good overlap |
| tune3, tunea2, dim | 5, 12, 5 | 500 / 500 / 200 | no script in repo; not in `keep_fam` of prep_cv_summary.R | superseded exploration |

Summary tables in `output/cv/basic_dgp/summaries/` (`cv_summary.csv` is stale: built from the
2 Sep tune*.rds). Figures via `analysis/cv/basic_dgp/plot_grids.R`.

### cv / wz_replication (`runs/cv/wz_replication/`, `results/cv/wz_replication/`, n = 5000, 1000 reps, master seed 20260903)

| batch | cells | script | grid |
|---|---|---|---|
| wc_overlap_c1..c4 | 4 | run_wc_overlap.R | logit scale c 1..4 |
| wc_noise_s{1,10,30,100}_c{1..4} | 16 | run_wc_noise.R | noise SD x overlap |
| wc_on_n{1000..20000}_c{1..4} | 20 | run_wc_overlap_n.R | n x overlap |
| wc_basis_K{10,20,65,125}_c{1,3} | 8 | run_wc_basis.R | basis size x overlap |
| wc_att_v2 | 1 | producer not in repo (run_wc_path.R writes wc_att_v1, absent) | ATT path, K = 20 |
| wc_sbw_v1, wc_sbw_v2 (+ _meta) | 2 | producer not in repo | sbw replication of Table 4 |
| wc_tuner_n{1000,5000,20000}_c{1..4} | 12 | producer not in repo; not read by any script | tuner batches |

## 4. Findings (verbatim)

### ipw: correct specification (`analysis/ipw/correctspec.qmd`)

> **ALL CV.BALNET ESTIMATES AND DIAGNOSTICS ARE INVALID UNTIL FIX.** (line 79)

> **Balnet (path endpoint) RMSE dominates all estimators including both oracles (HT and Hájek) at every** $\bf n$**,** consistent with Wager (2024). *(single exception: quad1 n=500 vs bal05 is unresolved, paired z = −0.13, a tie by our convention).* (line 83)

> **at** $\bf n=500$**,** balancing methods outperform MLE. (line 86)

> **as** $\bf n$ **increases,** penalized MLE `glmnet.cv` consistency property wins out. (line 87)

> `bal05` outperforms `glmnet` at $n \leq 1000$; `bal10` only at $n = 500$**.** As $n$ increases, `bal05` and `bal10` are outperformed by `glmnet` as their bias converges to a floor as their imbalance converges to their set maximum imbalance, $\lambda$. (line 88)

> **balnet.cv over-regularizes, erasing the loss function advantage.** (line 93)

> median $\lambda_{cv}$ decays 0.011 → 0.006 over $n$ = 500 → 50,000 ($\approx n^{-0.12}$), against $n^{-0.53}$ for `cv.glmnet` and $\lambda_{end}$ = 1e-4 throughout. In KKT terms, the selected maximum allowable imbalance barely shrinks with $n$. (line 96)

> Bal0's performance is invariant to sign changes, confirming that exact balancing is insensitive to confounding alignment (the degree to which propensity coefficient signs match outcome coefficient signs). (line 101, first sentence)

> Under correct specification non-normalized MLE *(glmHT)* and normalized MLE *(glmHajek)* are competitive. However, glmHajek imposes the Hájek normalization by construction (Chattopadhyay et al., 2020), and introducing mispeccification causes *glmHT's* RMSE and bias to explode, so we will only use the normalised version for our analysis. (line 117)

> Cross-validation targeting CB loss rejects the best-performing λ at small n, and is indifferent to it at large n. (line 197)

> 5 curves, one for each CV fold. Should be 1 CV curve. Drops occur every 20 indicies, for each CV fold. (line 244)

> confirmed our fix, taking the average over rows, produces 1 smooth curved row. (line 318)

### ipw: the cv.balnet fix (`analysis/ipw/cvfix.qmd`)

> Only `balnetcv` and its diagnostics (`lam_balcv`, `logcvloss_cv`, `logcvloss_end`, `smd_cv`, `nnz_balcv`) differ between builds; every other row is identical to full precision in both specs *(verified byte-identical on the raw files)*. Seed sharing worked, so any performance change below is attributable to the CV fix alone. (line 66)

> **the fix hurts at** $\bf n=500$ *(RMSE up \~25 to 35% in both specs)*, is roughly even at 1000, then wins decisively from 5000 on. fix/bug RMSE ratio at 50k: \~0.16 to 0.18 correct spec, \~0.47 to 0.71 misspec. (line 117)

> **the CV-rate critique was substantially a bug artifact.** Median $\lambda_{cv}$ decay under the fix is $\approx n^{-0.53}$ (correct) and $\approx n^{-0.50}$ (misspec), against $\approx n^{-0.11}$ on the bugged runs; it now matches `cv.glmnet`'s $\approx n^{-0.52}$. The E8 / CV-critique text (rate channel) needs rewriting before any of it goes in the draft. (line 201)

> still to check on the fixed runs before reusing the old claims: the **level** channel (is $\lambda_{cv}$ still above the RMSE-optimal $\lambda$, cite Wyss 2026 + Ertefaie for that channel only) and the **selection-noise variance** channel (does cor($\lambda_{cv}, \hat\tau$) still approach $-1$, and does bias$^2$ share still plateau below 100%). (line 202)

### ipw: overlap and misspecification (`analysis/ipw/overlap.qmd`)

> **summary:** our expectation that under weak overlap balancing estimators performance would degrade at a reduced rate compared to MLE holds. all n. (line 49)

> **Bal0 is the best estimator, consistently outperforming or matching other estimators.** (line 55)

> **Covariate balancing outperforms MLE over every overlap level, and for every** $n$**.** (line 63)

> **balnet.cv becomes less competitive against bal0 and bal05 as overlap degrades:** this calls into question its CV criterion. (line 93)

> **Path truncation mechanism: Bal05 is the only estimator competitive with Bal0, due to path truncation.** As overlap degradrs, the $\lambda_{endpoint}$ increases, and with increases in $n$ it decreases. Bal0 generally either outperforms or matches bal05, but in smaller samples as overlap degrades bal05 becomes increasingly competitive. (line 97)

> As expected, misspecification and overlap compound, and MLE methods are considerably outperformed by balancing methods. MLE RMSE increases substantially, with its bias behaving erratically as it converges to the wrong prediction (sign flips under quadratic), and extreme behaviour concentrated under quadratic outcomes, while the balancing estimators remain stable. (line 143)

> Within balnet, $\lambda$ is a bias/variance dial where less is better at every $n$ and overlap level; balcv provides no evidence for cross validation, even at small $n$ and weak overlap, where it over-regularizes most (quarantined pending the cv.balnet fix). In large samples, its bias floor is too high. Despite marginal wins\*, Still no strong evidence supporting for setting $\lambda > 0$. (line 147)

> **Extremes:** MLE produces 28 extreme estimates ($|\hat{\tau}| > 10$), concentrated in quadratic outcomes and moderate overlap ($0.5 \leq c \leq 1$). Balancing produces 0. *No need for ad hoc weight trimming/truncation with balnet.* (line 181)

> **balcv is third within balnet at** $c < 2$, and second at $c = 2$, where it beats bal05 at every $n$; it never beats bal0. (line 209)

### ipw: dimensionality (`analysis/ipw/Dimensionality.qmd`)

> As expected, balancing estimators are more robust to high-dimensionality thn MLE estimation and Bal0 is still the best performing estimator, converging with bal05's path endpoint as $p/n$ increases. The only exception is for balnet.cv, which eventually underperforms MLE as dimensionality increases. (line 55)

> These results support that the balnet.cv cross validation criterion *(out-of-sample covariate balancing loss)* is suboptimal. (line 61)

> balcv degrades *faster* than glmnet as $p/n$ grows and crosses below glm_hj . The loss function survives dimensionality, but the CV rule doesn't. **strong evidence against balnet.cv.** (line 97)

### ipw: sparsity, SNR, treated share (`increased_sparsity.qmd`, `SNR.qmd`, `unequal_treatment.qmd`)

> BAL0 \> MLE \> BALCV: (increased_sparsity.qmd line 81)

> bal0 by n = 10k, and active covariates = p = 50, reducing decay improves performance for MLE and balCV, they perform almost on par with bal0. bal05 and bal1q0 ared significantly degraded. (increased_sparsity.qmd line 87)

> bal0 outperforms again (SNR.qmd line 52)

> bal05, bal10 are outperformed by MLE with increases in SNR. (SNR.qmd line 56)

> Under quadratic outcomes: bal0 also degrades with SNR increases. (SNR.qmd line 60)

> Bal0 is the performs the best when treatment allocation is degraded. (unequal_treatment.qmd line 52)

> balnet.cv is outperformed by MLE, cross validation criterion is sensitive to unequal treatment allocations. (unequal_treatment.qmd line 54)

`Instruments.qmd`, `sparsity_hd.qmd`, `sparsity_hd_check.qmd` and `cv_analysis.qmd` contain no written findings yet.

### cv / basic_dgp (`analysis/cv/basic_dgp/cv_summary.qmd`)

> **A:** Is the optimal (minimizes RMSE) $\lambda$ interior, i.e. not the path floor $\lambda_{\min}$?\ (line 36)

> **B:** If so, does any selector beat the path floor? (line 37)

> **A:** with `sigma_y = 5`, interior solutions exist at all 3 overlap levels. (line 69)

> **B:** With large outcome noise, selectors outperform the path floor $\lambda_{\min}$ at good and moderate overlap; at bad overlap the best selector only ties it. (line 71)

> **Increasing outcome noise causes both interior optimums and selector wins.** (line 79)

> **Bootstrapping vs cross-validation:** for both criteria where the two are available (`smd`, `inf`), bootstrapping outperforms CV at every overlap level. `cv.bloss` is the only competitive criterion under bad overlap, and has no bootstrapped equivalent. (line 83)

> Under both high and low noise settings *(Appendix 3)* outcome density does not appear to drive selection performance, or interior minimums. Results for density show only minor changes to relative performance, with any ranking differences likely sampling noise. (line 99)

> **A:** interior minimums appear at lower noise under strong overlap; as overlap degrades, more noise is needed. (line 107)

> **Criterion selection:** `cv.bloss` and `boot.inf` are the most effective criteria. `cv.bloss` is more robust to overlap degradation, outperforming `boot.inf` under "bad" overlap for $\sigma_y \le 5$ and underperforming for "good" and "moderate" overlap. (line 129)

> **overlap axis** ($\sigma_y = 1$): there is no interior solution, the path is still falling, indicating $\lambda_{\min}$ is best at every level and its apparent losses in the figure (`c_prop` 6 to 16) are likely artifacts . The order of `cv.bloss` and `boot.inf` flips between `c_prop` 2 and 2.5. (line 145)

> **noise axis** (good overlap): $\lambda_{\min}$ is best only at $\sigma_y$ = 1; `boot.inf` ties it at 2 and beats it from 3, `cv.bloss` ties it at 5 and beats it from 7; the order of `cv.bloss` and `boot.inf` flips between $\sigma_y$ 7 and 10. (line 151)

> As selection models are outcome blind, increasing noise doesnt change their selection lambda, it only raises optimal level of imbalance, CV bloss selects higher, so when noise is high enough - it wins. (line 202)

> **Flip dynamic related to criterion not bootstrapping.** When outcome isnt visible, the `inf` (max $|\mathrm{SMD}_j|$) criterion variants, which minimise the worst balanced co-variate, always choose a higher lambda, making them more robust to highly degraded overlap, since as overlap weakens the optimal RMSE rises. (line 254)

> **A:** interior minimums appear at lower noise as $\alpha$ rises: only at $\sigma_y = 10$ for $\alpha = 0$, from $\sigma_y = 1$ for $\alpha = 0.25$ to 0.75. (line 346)

> **Bootstrapping vs cross-validation:** bootstrap \> CV everywhere except $\alpha = 0.75$, $\sigma_y = 10$, where `cv.inf` \> `boot.inf` and `boot.inf` drops below `boot.smd`. (line 352)

Earlier synthesis (`archive/exploring_cv/cv_summary.Rmd`, superseded draft):

> Against the lasso floor the verdict is: selectors win only where $R^2 \lesssim 0.1$ or overlap is very weak, by 4–8% at $R^2 = 0.08$ and 13–31% at $R^2 = 0.02$, and that advantage grows with $n$ under weak overlap. At $R^2 \ge 0.36$ the lasso floor wins by 15–57%, by a margin that grows with $n$, with $p$ and with confounder spread. At $p \ge 500$ the lasso floor beats every elastic-net selector at $R^2 \ge 0.2$; at $R^2 = 0.02$ selectors win by 7–23% at $p = 500$ and the gap closes at $p = 1000$. (line 72)

> All five selectors are outcome-blind, so the $\lambda$ they choose does not move with $\sigma_y$ while the optimum does. Each selector is a fixed rung; a selector "wins" a cell when the optimum has moved onto its rung. boot.inf sits nearest the optimum for $R^2 \ge 0.08$, cv.bloss for $R^2 \approx 0.02$. The mean-SMD criteria choose $\lambda_{\max}$ (no weighting) whenever $p \ge 500$ and are never competitive. (line 70)

### cv / wz_replication (`analysis/cv/wz_replication/wc_findings.qmd`)

> Replication: model B and both exact arms match Table 4 within 0.35; model Aapproximate does not (2.4 to 2.7 times theirs). Our loop picks a median\$\\delta = 0.02\$; theirs must pick tighter, or the code differs. Supplementcode not yet checked. (line 68)

> **Model A:** fails. RMSE more than doubles (0.92 to 2.11), all bias. Paper: 10% worse. (line 72)

> **Model B:** holds. RMSE falls 20% (0.29 to 0.23) - low SNR, perfect balancing fails. (line 74)

> The rule picks the same λ for both models; only the outcome, which it never sees, decides whether that λ helps. (line 76)

> **Model A: (equivalent to high SNR)** (line 158) Both claims fail. Best λ is the floor (0.43). (line 160) Every rule relaxes and loses; theirs by 3.6× (1.54). (line 162)

> Claim (1) holds: Best λ is 0.28 (0.19 oracle vs 0.32 at the floor). (line 166)

> Claim (2) fails: every rule stops far short of the oracle (λ = 0.28, RMSE 0.19). (line 168)

> The ranking is just the λ ranking, loser wins because B rewards not adjusting, on A it reverses. No rule can be right on both. (line 178)

## 5. Not yet run, failed, or missing

- `dimensionDGP2N5K` (`runs/ipw/DimensionDGP2N5K.R`, n 500/1000/5000 x p/n): only 5 of 12 cell
  checkpoints exist in `results/ipw/dimensionDGP2N5K_cells/`; no csv.gz; never in the ledger.
- `correctspecDGP200k` (first batch of `correctspec_dgp2_n800k.R`): no file; its output went to
  `correctspecDGP2.csv.gz` (old ledger row 31), which the 5-cell script cannot reproduce.
- `sparsityhighpDGP2` (old ledger rows 40-44): failed five times, superseded by `sparsityE6DGP2`.
- ipw batches whose scripts are not in the repo: `correctspecDGP2_ks / _pos / _mixed`, and the
  `cvtune*` batches (`results2/`) that `archive/analysis/cv_tuning.qmd` expected; never produced.
- cv / thesis_dgp: reserved in the layout, no scripts and no results.
- cv / basic_dgp families with results but no run script: `tunen2`, `dimsnr`, `snr_good` (used by
  the summaries), `tune3`, `tunea2`, `dim` (superseded); `tunep2` ran only at 50 reps (planned 200).
- cv / wz_replication: producers of `wc_att_v2`, `wc_sbw_v*`, `wc_tuner_*` are not in the repo;
  `run_wc_path.R` writes `wc_att_v1`, which does not exist; the sbw Supplement code check
  (line 68 above) is outstanding.
- Analysis still to write: `Instruments.qmd`, `sparsity_hd*.qmd`, `cv_analysis.qmd` have no
  prose; cvfix.qmd lists the level and selection-noise channels to re-check on the fixed runs
  and the misspecification effect on lambda_cv.
- Open items from the phase 4 report: 29 uncaptured warnings in the regression run; 4 empty
  outcome entries per column at n = 500 in both `correctspecDGP2` files.
- The stored `cv_summary.csv` predates the 10 Sep `tune*` reruns; `prep_cv_summary.R` and
  `plot_grids.R` have not been rerun on the restructured tree.
- `registry.csv` is empty since Phase 5; every batch above predates the new ledger.
