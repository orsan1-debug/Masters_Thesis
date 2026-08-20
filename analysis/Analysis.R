# ============================================================================
# analysis_setup.R -- shared setup for the E1-E6 analysis scripts.
# ============================================================================

library(dplyr)
library(tidyr)      
library(ggplot2)
library(patchwork)  # (l | r) panels
library(gt)
# pin dplyr verbs against masking (MASS, stats, ...)
select <- dplyr::select
filter <- dplyr::filter

drop0 <- function(x) sub("^(-?)0\\.", "\\1.", x)   # ".05" not "0.05" on axes

# ---- 1. Load ---------------------------------------------------------------

design_cols <- c("n", "p", "s", "decay", "misspec", "covcor", "outcome", "overlap",
                 "signs", "strength", "decay_ps", "decay_out", "outcome_set",
                 "treat_prop")
outcome_cols <- c("linear", "quad1", "quad2", "exp")

load_sim <- function(path) {
  d <- tibble::as_tibble(data.table::fread(path))
  d <- mutate(d, across(any_of(outcome_cols), \(x) suppressWarnings(as.numeric(x))))
  if ("err" %in% names(d)) d$err <- coalesce(as.character(d$err), "")
  pivot_longer(d, any_of(outcome_cols), names_to = "outcome", values_to = "tau_hat")
}

# diagnostics are logged as rows of `estimator`; identify them by name
is_diag <- function(x)
  grepl("^(lam_|nnz_|smd[01]_|cvloss_|logcvloss_|trunc05|prev|emin|emax|nout)", x) |
  x %in% c("trunc05", "prev", "emin", "emax", "nout05", "nout01")


# ---- 2. Cell summaries -----------------------------------------------------

# Screen n_na / n_extreme before reading anything else; bias/empse/rmse
# propagate NA loudly by design.
summarise_sim <- function(data, by = c(design_cols, "estimator")) {
  by <- intersect(by, names(data))
  data |>
    filter(!is_diag(estimator)) |>
    group_by(across(all_of(by))) |>
    summarise(
      m          = dplyr::n(),
      bias       = mean(tau_hat),
      empse      = sd(tau_hat),
      rmse       = sqrt(mean(tau_hat^2)),
      mcse_bias  = empse / sqrt(m),
      mcse_empse = empse / sqrt(2 * (m - 1)),
      mcse_rmse  = sd(tau_hat^2) / (2 * rmse * sqrt(m)),
      n_extreme  = sum(abs(tau_hat) > 10, na.rm = TRUE),
      n_na       = sum(is.na(tau_hat)),
      .groups    = "drop"
    )
}

# ---- 3. Estimator labels, palette, sets ------------------------------------

lab_est <- c(balnet0        = "BalNet (\u03bb=0)",
             balnet05       = "BalNet (\u03bb=.05)",
             balnet10       = "BalNet (\u03bb=.10)",
             balnetcv       = "BalNet (CV)",
             balnetrate     = "BalNet (\u221a(log p/n))",
             glmnetcv_hajek = "GLM (Hajek)",
             glmnetcv_ht    = "GLM (HT)",
             oracle_hajek   = "Norm. Oracle",
             oracle_ht      = "Oracle (HT)")
pal_est <- c(balnet0 = "#009E73", balnet05 = "#56B4E9", balnet10 = "#E69F00",
             balnetcv = "#0072B2", balnetrate = "#D55E00",
             glmnetcv_hajek = "#B22222", glmnetcv_ht = "#E7298A",
             oracle_hajek = "grey55", oracle_ht = "#404040")

est_main <- c("balnet0", "balnet05", "balnet10", "balnetcv", "glmnetcv_hajek")
est_bal  <- c("balnet0", "balnet05", "balnet10", "balnetcv")

# ---- 4. Theme, scales, facet labels ----------------------------------------

theme_sim <- theme_minimal(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.title     = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1))

scale_x_n <- scale_x_log10(breaks = c(500, 1000, 5000, 10000, 50000),
                           labels = c("500", "1k", "5k", "10k", "50k"))
lab_n   <- as_labeller(\(x) paste0("n = ", format(as.numeric(x), big.mark = ",")))
n_lab   <- c("1000" = "1k", "5000" = "5k", "10000" = "10k", "50000" = "50k")
oc_labs <- c(linear = "Linear", quad1 = "Quadratic", exp = "Exponential")
as_c <- function(x) factor(x, sort(unique(x), decreasing = TRUE), paste0("c = ", drop0(sort(unique(x), decreasing = TRUE))))

# shorten table headers ("5000" -> "5k"); anything unlisted passes through
n_relabel <- function(v) unname(ifelse(v %in% names(n_lab), n_lab[v], v))

# ---- 5. Core line plot + twin panel ----------------------------------------

# metric vs x, coloured by estimator; caller appends x-scale and facets
line_plot <- function(data, y, ylab, x = n, ests = est_main, from_zero = TRUE) {
  ests <- intersect(names(lab_est), ests)   # canonical order, no ghost keys
  p <- data |>
    filter(estimator %in% ests) |>
    ggplot(aes({{ x }}, {{ y }}, colour = estimator, group = estimator)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
    scale_colour_manual(values = pal_est, breaks = ests,
                        labels = \(b) unname(lab_est[b])) +
    labs(y = ylab, colour = NULL) + theme_sim
  if (from_zero) p + expand_limits(y = 0) else p
}

# RMSE | |Bias| side by side, one shared legend
panel <- function(data, oc, x = n, ests = est_main, extra = NULL,
                  title = oc_labs[oc]) {
  d <- filter(data, outcome == oc)
  l <- line_plot(d, rmse,      "RMSE",   x = {{ x }}, ests = ests) + extra
  r <- line_plot(d, abs(bias), "|Bias|", x = {{ x }}, ests = ests) + extra
  (l | r) + plot_layout(guides = "collect") +
    plot_annotation(title = title) & theme(legend.position = "bottom")
}

# ---- 6. MSE decomposition (bias^2 / variance shares) -----------------------

pal_comp <- c("Variance" = "#377eb8", "Bias\u00b2" = "#e41a1c")

decomp_plot <- function(data, ests = est_main, x = n, cols = outcome,
                        x_scale = scale_x_n, xlab = "n",
                        col_lab = labeller(outcome = oc_labs, n = lab_n),
                        title = NULL) {
  keep <- intersect(names(lab_est), ests)
  data |>
    filter(estimator %in% keep) |>
    mutate(denom    = bias^2 + empse^2,
           bias_sq  = bias^2  / denom,
           variance = empse^2 / denom,
           estimator = factor(lab_est[estimator], levels = lab_est[keep])) |>
    pivot_longer(c(variance, bias_sq), names_to = "component", values_to = "share") |>
    mutate(component = factor(component, c("variance", "bias_sq"),
                              c("Variance", "Bias\u00b2"))) |>
    ggplot(aes({{ x }}, share, fill = component)) +
    geom_area() +
    facet_grid(rows = vars(estimator), cols = vars({{ cols }}), switch = "y",
               labeller = col_lab) +
    x_scale +
    scale_y_continuous(breaks = 1, labels = scales::percent) +
    scale_fill_manual(values = pal_comp) +
    labs(x = xlab, y = "Share of MSE", fill = NULL, title = title) +
    theme_sim +
    theme(strip.text.y.left = element_text(angle = 0))
}

# ---- 7. Performance table --------------------------------------------------

# estimators x design columns, ordered by mean metric.
# cols = "n" -> flat; cols = c("overlap","n") -> c-spanners;
# group = "outcome" -> stacked row groups.
perf_table <- function(data, ests = est_main, cols = "n",
                       top = "rmse", bottom = "abs_bias", digits = 2,
                       order_by = "rmse", title = "", subtitle = NULL,
                       group = NULL,
                       span_lab = \(s) paste0("c = ", s), span_desc = TRUE) {
  stopifnot(length(cols) %in% 1:2,
            c(top, bottom, order_by) %in% c("rmse", "bias", "abs_bias", "empse"))
  keep <- intersect(names(lab_est), ests)
  cvar <- cols[length(cols)]
  svar <- if (length(cols) == 2L) cols[[1]] else NULL
  pm   <- c(rmse = "RMSE", bias = "Bias", abs_bias = "|Bias|", empse = "EmpSE")
  f    <- paste0("%.", digits, "f")
  grp  <- !is.null(group)
  
  d <- data |> filter(estimator %in% keep) |> mutate(abs_bias = abs(bias)) |> mutate(
    Estimator = factor(unname(lab_est[estimator]), levels = unname(lab_est[keep])),
    cell = sprintf(paste0(f, "<br><small>(", f, ")</small>"), .data[[top]], .data[[bottom]]),
    avg  = ave(.data[[order_by]], estimator, FUN = mean),
    .col = if (is.null(svar)) paste0(cvar, "=", .data[[cvar]])
    else paste0(svar, "=", .data[[svar]], "|", cvar, "=", .data[[cvar]]))
  if (grp) d <- mutate(d, Group = factor(oc_labs[.data[[group]]], levels = oc_labs))
  
  cv  <- sort(unique(d[[cvar]]))
  sv  <- if (is.null(svar)) NULL else sort(unique(d[[svar]]), decreasing = span_desc)
  key <- if (is.null(svar)) paste0(cvar, "=", cv)
  else c(t(outer(sv, cv, \(s, c) paste0(svar, "=", s, "|", cvar, "=", c))))
  
  id_vars <- c(if (grp) "Group", "Estimator", ".col", "cell", "avg")
  wide <- d |> distinct(across(all_of(id_vars))) |>
    pivot_wider(names_from = .col, values_from = cell)
  wide <- if (grp) arrange(wide, Group, avg) else arrange(wide, avg)
  present <- intersect(key, names(wide))
  wide <- wide |> select(any_of("Group"), Estimator, all_of(present))
  
  inner_lab <- \(x) n_relabel(sub("^[^=]*=", "", sub(".*\\|", "", x)))
  if (is.null(subtitle)) subtitle <- paste0(pm[[top]], " (", pm[[bottom]], ")")
  
  g <- if (grp) gt(wide, groupname_col = "Group") else gt(wide)
  g <- g |> tab_header(title = title, subtitle = md(subtitle))
  for (s in sv) g <- tab_spanner(g, span_lab(s),
                                 starts_with(paste0(svar, "=", s, "|")))
  g |>
    cols_label_with(all_of(present), inner_lab) |>
    fmt_markdown(all_of(present)) |> cols_align("center", all_of(present)) |>
    tab_style(cell_text(weight = "bold"), cells_body(Estimator)) |>
    opt_table_lines("default") |>
    tab_source_note(md(sprintf("_Top %s / bottom %s; sorted by mean %s._",
                               pm[[top]], pm[[bottom]], pm[[order_by]])))
}

# ---- 8. Path-truncation table ----------------------------------------------

# % reps where the attained path endpoint sits above lambda = .05 in either
# arm (trunc05 as logged) -- the positivity diagnostic.
trunc_table <- function(data, by = c("overlap", "n"), cols = by[length(by)],
                        outcome_name = NULL, title = NULL) {
  dw <- as_diag(data)
  by <- intersect(by, names(dw))
  stopifnot(length(cols) == 1, cols %in% by)
  wide <- dw |>
    filter(if (is.null(outcome_name)) TRUE else outcome == outcome_name) |>
    group_by(across(all_of(by))) |>
    summarise(pct = sprintf("%.1f%%", 100 * mean(trunc05, na.rm = TRUE)),
              .groups = "drop") |>
    pivot_wider(names_from = all_of(cols), values_from = pct,
                names_prefix = paste0(cols, "="))
  ord  <- grep(paste0("^", cols, "="), names(wide), value = TRUE)
  ord  <- ord[order(as.numeric(sub(".*=", "", ord)))]
  wide <- select(wide, all_of(setdiff(names(wide), ord)), all_of(ord))
  if ("overlap" %in% names(wide)) names(wide)[names(wide) == "overlap"] <- "c"
  g <- gt(wide)
  if (!is.null(title)) g <- tab_header(g, title = md(title))
  g |>
    cols_label_with(all_of(ord), \(x) n_relabel(sub(".*=", "", x))) |>
    cols_align("center")
}

# ---- 9. Misspecification helpers -------------------------------------------

# delta = misspecified - correct at matching (n, outcome, overlap, estimator)
delta_ms <- function(summ, ests = est_main) {
  summ |>
    filter(estimator %in% ests) |>
    mutate(absbias = abs(bias)) |>
    select(n, outcome, overlap, estimator, misspec, rmse, absbias) |>
    pivot_wider(names_from = misspec, values_from = c(rmse, absbias)) |>
    mutate(d_rmse = rmse_TRUE - rmse_FALSE,
           d_bias = absbias_TRUE - absbias_FALSE)
}

delta_table <- function(d, ests = est_main, at_n = 5000, ov = c(0.25, 0.5, 0.75, 1),
                        digits = 3,
                        title = paste0("\u0394 (misspecified \u2212 correct), n = ",
                                       format(at_n, big.mark = ","))) {
  keep <- intersect(names(lab_est), ests)
  f    <- paste0("%.", digits, "f")
  wide <- d |>
    filter(n == at_n, overlap %in% ov, estimator %in% keep) |>
    mutate(outcome   = factor(oc_labs[outcome], oc_labs),
           Estimator = factor(unname(lab_est[estimator]), unname(lab_est[keep])),
           cell = sprintf(paste0(f, "<br><small>(", f, ")</small>"),
                          d_rmse, d_bias)) |>
    arrange(Estimator) |>
    select(outcome, overlap, Estimator, cell) |>
    pivot_wider(names_from = Estimator, values_from = cell) |>
    arrange(outcome, desc(overlap))
  est_cols <- unname(lab_est[keep])
  gt(wide, groupname_col = "outcome", rowname_col = "overlap") |>
    tab_header(title = md(title), subtitle = md("\u0394RMSE (\u0394|Bias|)")) |>
    tab_stubhead(label = "c") |>
    fmt_markdown(all_of(est_cols)) |>
    cols_align("center", all_of(est_cols)) |>
    tab_source_note(md("_Top \u0394RMSE / bottom \u0394|Bias|._"))
}


# deltaRMSE | delta|Bias| twin vs overlap, faceted by n, walked over outcomes
delta_panel <- function(data, oc, ests = est_main, ov = c(0.25, 0.5, 0.75, 1),
                        title = oc_labs[oc]) {
  d <- filter(data, outcome == oc, overlap %in% ov)
  shared <- list(
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4),
    scale_x_reverse(breaks = ov),
    facet_wrap(~n, labeller = lab_n),
    scale_y_continuous(labels = drop0))
  l <- line_plot(d, d_rmse, NULL, x = overlap, ests = ests, from_zero = FALSE) +
    ggtitle("\u0394RMSE") + shared
  r <- line_plot(d, d_bias, NULL, x = overlap, ests = ests, from_zero = FALSE) +
    ggtitle("\u0394|Bias|") + shared
  (l | r) + plot_layout(guides = "collect") +
    plot_annotation(title = title) & theme(legend.position = "bottom")
}

# |Bias| and SD vs overlap, coloured by specification
pal_spec <- c(Correct = "#0072B2", Misspecified = "#B22222")

spec_grid <- function(data, ests = est_bal, at_n = 5000, ov = c(0.25, 0.5, 0.75, 1, 2)) {
  ests <- intersect(names(lab_est), ests)
  data |>
    filter(estimator %in% ests, n == at_n, overlap %in% ov) |>
    pivot_longer(c(bias, empse), names_to = "metric") |>
    mutate(metric    = factor(metric, c("bias", "empse"), c("|Bias|", "SD")),
           spec      = ifelse(misspec, "Misspecified", "Correct"),
           estimator = factor(lab_est[estimator], levels = lab_est[ests]),
           outcome   = factor(outcome, names(oc_labs), oc_labs)) |>
    ggplot(aes(overlap, abs(value), colour = spec)) +
    geom_line() + geom_point(size = 1.5) +
    facet_grid(metric ~ estimator + outcome, scales = "free_y") +
    scale_x_continuous(breaks = c(0.25, 0.5, 1, 2), labels = c(".25", ".5", "1", "2")) +
    scale_colour_manual(values = pal_spec) +
    labs(x = "Overlap (c)", y = NULL, colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          strip.text.x = element_text(size = 7),
          axis.text.x  = element_text(angle = 45, hjust = 1, size = 7))
}

# ---- 10. Extreme-estimate heatmap ------------------------------------------

# reps with |tau_hat| > 10, overlap x n, faceted by outcome. Tile + numerator
# = the estimator's count per cell; denominator = its total across the figure.
extreme_heat <- function(data, est = "glmnetcv_hajek", x = overlap,
                         xlab = "Overlap (c)", x_desc = TRUE,
                         ov = NULL, title = NULL) {
  d <- filter(data, estimator == est, if (is.null(ov)) TRUE else {{ x }} %in% ov)
  total <- sum(d$n_extreme)
  if (is.null(title))
    title <- sprintf("%s: extreme estimates (|\u03c4\u0302| > 10)", unname(lab_est[est]))
  d |>
    group_by(outcome, n, x = {{ x }}) |>
    summarise(extreme = sum(n_extreme), .groups = "drop") |>
    mutate(fill_val = if_else(extreme > 0, extreme, NA_real_),
           outcome  = factor(outcome, names(oc_labs), oc_labs),
           n        = factor(n, sort(unique(n))),
           x        = factor(x, sort(unique(x), decreasing = x_desc))) |>
    ggplot(aes(x, n, fill = fill_val)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = if_else(extreme > 0, sprintf("%.0f/%.0f", extreme, total), "")),
              size = 3) +
    facet_wrap(~outcome) +
    scale_fill_gradient(low = "white", high = "#B22222", name = "Extreme\n(count)",
                        na.value = "white") +
    scale_y_discrete(labels = n_lab) +
    labs(x = xlab, y = "n", title = title) +
    theme_sim +
    theme(panel.grid = element_blank(), legend.position = "right")
}

# ---- 11. Diagnostics infrastructure ----------------------------------------

# Diagnostics arrive as estimator rows; diag_wide() lifts them to one row per
# rep (per outcome -- values repeat across outcomes). as_diag() accepts raw
# or already-wide input, so every helper below can be fed replicate data
# directly.
diag_wide <- function(data) {
  id <- intersect(c(design_cols, "sim"), names(data))
  d  <- filter(data, is_diag(estimator))
  if (!nrow(d)) stop("no diagnostic rows found")
  dw <- pivot_wider(d, id_cols = all_of(id), names_from = estimator,
                    values_from = tau_hat)
  err_tab <- data |> distinct(across(all_of(id)), err) |>
    mutate(err = as.integer(err != ""))
  left_join(dw, err_tab, by = id)
}

as_diag <- function(data)
  if ("lam_glmcv" %in% names(data)) data else diag_wide(data)

# med / q10 / q90 columns -> long, one row per stat x quantile set
pivot_q <- function(dg, prefix)
  pivot_longer(dg, matches(paste0("^", prefix, "_(med|q10|q90)$")),
               names_to = c("stat", ".value"),
               names_pattern = paste0("^(", prefix, ")_(med|q10|q90)$"))

# cell-level med / q10 / q90 of every logged quantity + failure counts
summarise_diag <- function(data, by = design_cols) {
  dw <- as_diag(data)
  by <- intersect(by, names(dw))
  dw |>
    group_by(across(all_of(by))) |>
    summarise(
      m         = dplyr::n(),
      n_err     = sum(err != 0, na.rm = TRUE),
      p_trunc05 = mean(trunc05, na.rm = TRUE),
      across(matches("^(lam|nnz|smd|cvloss|logcvloss|prev|emin|emax|nout)"),
             list(med = \(x) median(x, na.rm = TRUE),
                  q10 = \(x) quantile(x, .10, na.rm = TRUE, names = FALSE),
                  q90 = \(x) quantile(x, .90, na.rm = TRUE, names = FALSE))),
      .groups = "drop"
    )
}

# ---- 12. Lambda paths ------------------------------------------------------

# Median (optional q10-q90 ribbon) of CV-selected lambda, attained path
# endpoint, and glmnet lambda.min vs n. Dashed black = sqrt(log p / n)
# (Wager 2024, s.7.2); dotted = the fixed-lambda marks. Colour = quantity,
# linetype = arm (the CB loss fits one model per arm).
lam_meta <- tibble::tribble(
  ~stat,        ~quantity,           ~arm,
  "lam_balcv1", "CV \u03bb",         "treated",
  "lam_balcv0", "CV \u03bb",         "control",
  "lam_end1",   "path endpoint",     "treated",
  "lam_end0",   "path endpoint",     "control",
  "lam_glmcv",  "glmnet \u03bb.min", "single fit"
)
pal_lam <- c("CV \u03bb" = "#0072B2", "path endpoint" = "#009E73",
             "glmnet \u03bb.min" = "#B22222")
lt_arm  <- c(treated = "solid", control = "42", "single fit" = "solid")
lambda_plot <- function(dg, x = n, stats = lam_meta$stat,
                        ribbon = FALSE, rate = TRUE, ref_fixed = c(.05, .10),
                        p_dim = NULL) {
  # rate reference assumes x = n; set rate = FALSE for any other x
  d <- pivot_q(dg, "lam_[a-z0-9]+") |>
    filter(stat %in% stats) |>
    left_join(lam_meta, by = "stat")
  plt <- ggplot(d, aes({{ x }}, med, colour = quantity, linetype = arm,
                       group = interaction(quantity, arm)))
  if (ribbon)
    plt <- plt + geom_ribbon(aes(ymin = q10, ymax = q90, fill = quantity,
                                 group = interaction(quantity, arm)),
                             alpha = .12, colour = NA)
  if (rate) {
    rd <- dg |> distinct(across(any_of(design_cols)))
    if (!"p" %in% names(rd)) {
      if (is.null(p_dim))
        stop("rate line needs `p`: not in `dg`; supply p_dim or set rate = FALSE")
      rd$p <- p_dim
    }
    stopifnot(is.numeric(rd$p), is.numeric(rd$n))
    rd <- mutate(rd, med = sqrt(log(p) / n))
    plt <- plt + geom_line(data = rd, aes({{ x }}, med), inherit.aes = FALSE,
                           linetype = "dashed", colour = "black", linewidth = .5)
  }
  if (!is.null(ref_fixed))
    plt <- plt + geom_hline(yintercept = ref_fixed, linetype = "dotted",
                            colour = "grey45")
  plt +
    geom_line(linewidth = .7) + geom_point(size = 1.5) +
    scale_colour_manual(values = pal_lam) +
    scale_fill_manual(values = pal_lam, guide = "none") +
    scale_linetype_manual(values = lt_arm) +
    scale_y_log10() +
    labs(y = "\u03bb", colour = NULL, linetype = NULL) + theme_sim
}

# ---- 13. KKT / attained-imbalance check (balnet eq 10) ---------------------

# Bound: attained max|SMD| <= lambda per arm. Reports worst gap and share of
# reps exceeding it by > tol for whichever of {endpoint, .05, .10, CV} the
# file logged. Gaps <= 0 with slack shrinking in n are the fixed-lambda
# plateau mechanism. NB the driver logs a self-normalised SMD; small positive
# gaps may be that definitional mismatch, not solver failure.
kkt_check <- function(data, tol = 1e-6, by = design_cols) {
  dw  <- as_diag(data)
  by  <- intersect(by, names(dw))
  has <- \(...) all(c(...) %in% names(dw))
  if (has("smd1_0", "lam_end1"))
    dw$gap0  <- pmax(dw$smd1_0 - dw$lam_end1, dw$smd0_0 - dw$lam_end0)
  if (has("smd1_05")) dw$gap05 <- pmax(dw$smd1_05, dw$smd0_05) - 0.05
  if (has("smd1_10")) dw$gap10 <- pmax(dw$smd1_10, dw$smd0_10) - 0.10
  if (has("smd1_cv", "lam_balcv1"))
    dw$gapcv <- pmax(dw$smd1_cv - dw$lam_balcv1, dw$smd0_cv - dw$lam_balcv0)
  dw |>
    group_by(across(all_of(by))) |>
    summarise(m = dplyr::n(),
              across(starts_with("gap"),
                     list(max      = \(x) max(x, na.rm = TRUE),
                          pct_viol = \(x) 100 * mean(x > tol, na.rm = TRUE))),
              .groups = "drop")
}

# attained max|SMD| vs its ceiling: med (q10-q90) by arm; hline at the fixed
# lambda for at = "05"/"10". For at = "cv"/"0" the bound varies per rep --
# read jointly with kkt_check / lambda_plot.
smd_plot <- function(dg, x = n, at = "05", ribbon = TRUE) {
  d <- pivot_q(dg, paste0("smd[01]_", at)) |>
    mutate(arm = ifelse(startsWith(stat, "smd1"), "treated", "control"))
  p <- ggplot(d, aes({{ x }}, med, linetype = arm, group = arm))
  if (ribbon) p <- p + geom_ribbon(aes(ymin = q10, ymax = q90, group = arm),
                                   alpha = .12, colour = NA)
  if (at %in% c("05", "10"))
    p <- p + geom_hline(yintercept = as.numeric(paste0("0.", at)),
                        linetype = "dotted", colour = "grey45")
  p + geom_line(colour = "#0072B2", linewidth = .7) +
    geom_point(colour = "#0072B2", size = 1.5) +
    scale_linetype_manual(values = c(treated = "solid", control = "42")) +
    labs(y = paste0("attained max |SMD| at ",
                    switch(at, cv = "CV \u03bb", "0" = "path endpoint",
                           paste0("\u03bb = .", at))),
         linetype = NULL) + theme_sim
}

# ---- 14. Active-set size ---------------------------------------------------

lab_nnz <- c(nnz_balcv1 = "BalNet CV (treated)",
             nnz_balcv0 = "BalNet CV (control)",
             nnz_glm    = "glmnet \u03bb.min")

nnz_plot <- function(dg, x = n) {
  dg |>
    pivot_longer(matches("^nnz_[a-z0-9]+_med$"),
                 names_to = "stat", values_to = "med") |>
    mutate(stat = sub("_med$", "", stat)) |>
    ggplot(aes({{ x }}, med, colour = stat, group = stat)) +
    geom_line(linewidth = .7) + geom_point(size = 1.5) +
    scale_colour_manual(values = c(nnz_balcv1 = "#0072B2",
                                   nnz_balcv0 = "#56B4E9",
                                   nnz_glm    = "#B22222"),
                        labels = \(b) unname(lab_nnz[b])) +
    labs(y = "non-zero coefficients", colour = NULL) + theme_sim
}

# ---- 15. CV vs fixed lambda (E8: where and why CV loses) -------------------

# Per-rep join of selected estimates with that rep's lambdas. lam_cv averages
# the two arms because the bias mechanism sums both (the 2-lambda KKT term).
rep_cv <- function(data, oc, ests) {
  id <- intersect(c(setdiff(design_cols, "outcome"), "sim"), names(data))
  est <- data |>
    filter(estimator %in% ests, outcome == oc) |>
    select(all_of(id), estimator, tau_hat) |>
    pivot_wider(names_from = estimator, values_from = tau_hat)
  lam <- as_diag(data) |>
    filter(outcome == oc) |>
    transmute(across(all_of(id)),
              lam_cv  = (lam_balcv1 + lam_balcv0) / 2,
              lam_end = (lam_end1 + lam_end0) / 2)
  left_join(est, lam, by = id)
}

# WHERE: cell RMSE against the lambda each estimator sits at -- balnet0 at
# the median attained endpoint, .05/.10 fixed, balnetrate at sqrt(log p / n).
# balnetcv enters as a horizontal line, its selected-lambda q10-q90 mass as a
# band on the x axis: CV beats fixed lambda only if the line dips under the
# curve where its band sits. TODO once simulate_grid saves tau_path, replace
# these four points with the full tau(lambda) risk profile.
risk_curve <- function(data, oc, facets = vars(n), labeller = lab_n) {
  id <- intersect(setdiff(design_cols, "outcome"), names(data))
  q  <- rep_cv(data, oc, ests = "balnetcv") |>
    filter(!is.na(lam_cv)) |>                       # failed reps
    group_by(across(all_of(id))) |>
    summarise(cv_med  = median(lam_cv),
              cv_q10  = quantile(lam_cv, .10, names = FALSE),
              cv_q90  = quantile(lam_cv, .90, names = FALSE),
              end_med = median(lam_end), .groups = "drop")
  s <- summarise_sim(data) |> filter(outcome == oc) |> left_join(q, by = id)
  pts <- s |>
    filter(estimator %in% c("balnet0", "balnet05", "balnet10", "balnetrate")) |>
    mutate(lam = case_when(estimator == "balnet0"    ~ end_med,
                           estimator == "balnet05"   ~ 0.05,
                           estimator == "balnet10"   ~ 0.10,
                           estimator == "balnetrate" ~ sqrt(log(p) / n)))
  cv <- filter(s, estimator == "balnetcv")
  ggplot(pts, aes(lam, rmse)) +
    geom_rect(data = cv, aes(xmin = cv_q10, xmax = cv_q90, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = pal_est[["balnetcv"]], alpha = .10) +
    geom_vline(data = cv, aes(xintercept = cv_med),
               colour = pal_est[["balnetcv"]], linetype = "dashed") +
    geom_hline(data = cv, aes(yintercept = rmse), colour = pal_est[["balnetcv"]]) +
    geom_line(colour = "grey60") +
    geom_point(aes(colour = estimator), size = 2) +
    scale_colour_manual(values = pal_est, labels = \(b) unname(lab_est[b])) +
    scale_x_log10() +
    facet_wrap(facets, labeller = labeller, scales = "free_y") +
    labs(x = "\u03bb", y = "RMSE",
         caption = "line = BalNet (CV) RMSE; band = selected-\u03bb q10\u2013q90") +
    theme_sim
}

# WHY (level): performance of the CV estimate conditional on its selected
# lambda -- reps binned on lam_cv within cell, binned bias with an SD ribbon.
# Under linear / fixed4 / all-positive signs the large-n KKT ceiling gives
# bias ~ -2 sum(a_j) lambda = -5 lambda; pass ref_slope = -5 to overlay.
cv_bias_bins <- function(data, oc, bins = 10, ref_slope = NULL,
                         facets = vars(n), labeller = lab_n) {
  id <- intersect(setdiff(design_cols, "outcome"), names(data))
  d <- rep_cv(data, oc, ests = "balnetcv") |>
    filter(!is.na(lam_cv), !is.na(balnetcv)) |>
    group_by(across(all_of(id))) |>
    mutate(bin = dplyr::ntile(lam_cv, bins)) |>
    group_by(across(all_of(c(id, "bin")))) |>
    summarise(lam = median(lam_cv), bias = mean(balnetcv), sd = sd(balnetcv),
              .groups = "drop")
  p <- ggplot(d, aes(lam, bias)) +
    geom_ribbon(aes(ymin = bias - sd, ymax = bias + sd), alpha = .15) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = .4) +
    geom_line(colour = pal_est[["balnetcv"]]) +
    geom_point(colour = pal_est[["balnetcv"]], size = 1.5) +
    facet_wrap(facets, labeller = labeller, scales = "free") +
    labs(x = "selected \u03bb (bin median)", y = "mean \u03c4\u0302 within \u03bb-bin") +
    theme_sim
  if (!is.null(ref_slope))
    p <- p + geom_abline(intercept = 0, slope = ref_slope, linetype = "dotted")
  p
}

# WHY (dispersion): per-rep squared-error regret of CV against a fixed-lambda
# comparator, binned by lam_cv. True tau = 0, so tau_hat^2 is the per-rep
# squared error; regret > 0 means the CV pick did worse on that rep. Shows
# which part of the lam_cv distribution carries the loss.
cv_regret <- function(data, oc, ref = "balnet05", bins = 10,
                      facets = vars(n), labeller = lab_n) {
  id <- intersect(setdiff(design_cols, "outcome"), names(data))
  d <- rep_cv(data, oc, ests = c("balnetcv", ref)) |>
    filter(!is.na(lam_cv), !is.na(balnetcv)) |>
    mutate(regret = balnetcv^2 - .data[[ref]]^2) |>
    group_by(across(all_of(id))) |>
    mutate(bin = dplyr::ntile(lam_cv, bins)) |>
    group_by(across(all_of(c(id, "bin")))) |>
    summarise(lam = median(lam_cv), regret = mean(regret), .groups = "drop")
  ggplot(d, aes(lam, regret)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = .4) +
    geom_line(colour = "grey60") +
    geom_point(aes(colour = regret > 0), size = 1.8) +
    scale_colour_manual(values = c(`TRUE` = "#B22222", `FALSE` = "#009E73"),
                        labels = c(`TRUE` = "CV worse", `FALSE` = "CV better")) +
    facet_wrap(facets, labeller = labeller, scales = "free") +
    labs(x = "selected \u03bb (bin median)",
         y = paste0("mean \u03c4\u0302\u00b2 regret vs ", unname(lab_est[ref]))) +
    theme_sim
}

# ---- 16. CV criterion curves (what the criterion sees) ---------------------

# simulate_grid stores the full CV loss curve for the first cv_curve_reps
# reps of every cell in <out>_cvcurves.rds. Flat curves near the minimum are
# the candidate mechanism for the slow, dispersed lambda_cv (E8).
load_cv_curves <- function(out_file) {
  d      <- tibble::as_tibble(data.table::fread(out_file))
  cells  <- distinct(d[, seq_len(match("sim", names(d)) - 1L)])  # grid order
  curves <- readRDS(sub("\\.csv\\.gz$", "_cvcurves.rds", out_file))
  stopifnot(nrow(cells) == length(curves))
  bind_rows(lapply(seq_along(curves), function(ci) {
    cc <- bind_rows(lapply(seq_along(curves[[ci]]), function(r) {
      x <- curves[[ci]][[r]]
      if (is.null(x)) return(NULL)                  # failed / unlogged rep
      bind_rows(lapply(c("treated", "control"), function(a)
        tibble(rep = r, arm = a, lambda = x$lambda[[a]], cv = x$cv.mean[[a]])))
    }))
    if (!nrow(cc)) return(NULL)
    bind_cols(cells[rep(ci, nrow(cc)), ], cc)
  }))
}

# spaghetti of per-rep CV loss curves; black points mark each rep's selected
# lambda. Filter to one arm / cell subset before plotting, facet as needed.
cv_curve_plot <- function(cur, max_reps = 25, ref_fixed = c(.05, .10)) {
  cur <- filter(cur, rep <= max_reps)
  grp <- setdiff(names(cur), c("lambda", "cv"))
  sel <- cur |> group_by(across(all_of(grp))) |>
    slice_min(cv, n = 1, with_ties = FALSE) |> ungroup()
  ggplot(cur, aes(lambda, cv, group = interaction(rep, arm))) +
    geom_line(alpha = .2, colour = pal_est[["balnetcv"]]) +
    geom_point(data = sel, colour = "black", size = .7) +
    geom_vline(xintercept = ref_fixed, linetype = "dotted", colour = "grey45") +
    scale_x_log10() +
    labs(x = "\u03bb", y = "CV loss") + theme_sim
}

# ---- 17. Criterion insensitivity (E8: what the CV loss can't see) ----------

# Per-rep, per-arm shape of the CV loss surface from the sidecar curves.
# rel excess r(lambda) = cv(lambda)/min(cv) - 1: how much worse the criterion
# rates lambda than its own optimum. Flat region = {lambda : r < eps}; if the
# RMSE-best rung sits inside it, no selection rule on this loss can find it.
# Reference losses interpolated linearly on log-lambda; NA outside the grid.
curve_stats <- function(cur, ref_fixed = c(.05, .10), eps = .05) {
  grp <- setdiff(names(cur), c("lambda", "cv"))
  cur |>
    filter(is.finite(cv)) |>
    group_by(across(all_of(grp))) |>
    group_modify(\(d, key) {
      fmin <- min(d$cv)
      r_at <- \(l) approx(log(d$lambda), d$cv, xout = log(l), rule = 1)$y / fmin - 1
      flat <- d$lambda[d$cv / fmin - 1 < eps]
      tibble(lam_min   = d$lambda[which.min(d$cv)],
             r_end     = d$cv[which.min(d$lambda)] / fmin - 1,  # attained endpoint
             r_05      = r_at(ref_fixed[1]),
             r_10      = r_at(ref_fixed[2]),
             flat_lo   = min(flat), flat_hi = max(flat),
             flat_span = log10(max(flat) / min(flat)))          # decades wide
    }) |> ungroup()
}

# Cell summary: how flat, and which reference lambdas the criterion can't
# distinguish from its optimum. pct_* = share of stored reps with r < eps.
flatness_table <- function(cs, eps = .05, by = c(intersect(design_cols, names(cs)), "arm")) {
  cs |>
    group_by(across(all_of(by))) |>
    summarise(m_stored  = dplyr::n(),
              med_span  = median(flat_span),
              med_r_end = median(r_end),
              med_r_05  = median(r_05, na.rm = TRUE),
              med_r_10  = median(r_10, na.rm = TRUE),
              pct_end   = 100 * mean(r_end < eps),
              pct_05    = 100 * mean(r_05 < eps, na.rm = TRUE),
              pct_10    = 100 * mean(r_10 < eps, na.rm = TRUE),
              .groups = "drop")
}

# All-rep endpoint-vs-selected criterion gap from the logged two-point loss.
# Endpoint fits degenerate at small n (loss O(1e40+)), so log10 and tail
# shares only; the blowup IS the small-n endpoint-divergence evidence.
loss_gap <- function(data, by = design_cols) {
  as_diag(data) |>
    mutate(lr = pmax(log10(cvloss_end1 / cvloss_cv1),
                     log10(cvloss_end0 / cvloss_cv0))) |>
    group_by(across(all_of(intersect(by, names(data))))) |>
    summarise(m        = dplyr::n(),
              n_na     = sum(is.na(lr)),        # overflowed endpoint fits
              med_lr   = median(lr, na.rm = TRUE),
              q90_lr   = quantile(lr, .90, names = FALSE, na.rm = TRUE),
              pct_10x  = 100 * mean(lr > 1, na.rm = TRUE),
              pct_1e3x = 100 * mean(lr > 3, na.rm = TRUE),
              .groups = "drop")
}
