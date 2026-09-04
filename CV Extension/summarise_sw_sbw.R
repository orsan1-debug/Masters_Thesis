# summarise_wc_sbw.R -------------------------------------------------------
# Table 4 "Variance" rows (Wang & Zubizarreta 2020) from
# results/<batch_id>.csv.gz, one summary per estimand x truth x outcome x
# balance. Analysis only; reads results/.
# Population ATT for model A by Stein's lemma: eta = -Z1 - 0.1 Z4 is
# N(0, 1.01), E[s(eta)] = 1/2, and E[Zj s(eta)] = Cov(Zj, eta) E[s'(eta)],
# so E[Z1 | W = 1] = -2c and E[Z4 | W = 1] = -0.2c with
# c = E[s(eta){1 - s(eta)}]. Then PATT = 1.5(27.4(-2c) + 13.7(-0.2c))
# = -86.31c, about -17.8. PATE = 0. Model B has no effect.
# MCSE for RMSE uses the distribution-free MSE formula (Morris et al. 2019,
# Table 6), needed because model B is heavy-tailed.

library(dplyr)

dir      <- "C:/Users/otisr/Documents/Thesis 2026/Masters_Thesis/CV Extension/"
batch_id <- "wc_sbw_v2"
res <- data.table::fread(file = paste0(dir, "results/", batch_id, ".csv.gz"))

c0   <- integrate(\(x) plogis(x) * (1 - plogis(x)) * dnorm(x, 0, sqrt(1.01)),
                  -Inf, Inf)$value
patt <- -86.31 * c0

tab <- res |>
  mutate(
    population = case_when(outcome == "B" ~ 0, estimand == "ate" ~ 0,
                           TRUE ~ patt),
    sample     = case_when(outcome == "B" ~ 0, estimand == "ate" ~ sate,
                           TRUE ~ satt)
  ) |>
  tidyr::pivot_longer(c(population, sample), names_to = "truth",
                      values_to = "tau") |>
  mutate(err = tau_hat - tau) |>
  group_by(estimand, truth, outcome, balance) |>
  summarise(
    n_sim     = n(),
    bias      = mean(err),
    sd        = sd(err),
    rmse      = sqrt(mean(err^2)),
    rmse_mcse = sd(err^2) / sqrt(n()) / (2 * sqrt(mean(err^2))),
    delta     = mean(delta),
    .groups   = "drop"
  ) |>
  arrange(estimand, truth, outcome, balance)

print(tab, n = Inf)
table(res$delta[res$balance == "approx"], res$estimand[res$balance == "approx"])