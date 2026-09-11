# ============================================================================
# R/plots.R -- figure and table helpers, split verbatim from analysis/Analysis.R
# (Phase 3). Requires R/summarise.R to be sourced first: libraries, drop0(),
# design_cols, summarise_sim(), as_diag(), pivot_q(), rep_cv().
# ============================================================================

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
#' Overlap values as a descending "c = " factor
#'
#' @param x Numeric overlap values.
#' @return A factor with levels in decreasing order, labelled "c = <value>".
as_c <- function(x) factor(x, sort(unique(x), decreasing = TRUE), paste0("c = ", drop0(sort(unique(x), decreasing = TRUE))))

# shorten table headers ("5000" -> "5k"); anything unlisted passes through
#' Shorten sample-size labels for table headers
#'
#' @param v Character vector of values.
#' @return v with 1000 / 5000 / 10000 / 50000 replaced by 1k / 5k / 10k / 50k;
#'   other values pass through.
n_relabel <- function(v) unname(ifelse(v %in% names(n_lab), n_lab[v], v))

# ---- 5. Core line plot + twin panel ----------------------------------------

# metric vs x, coloured by estimator; caller appends x-scale and facets
#' Metric against a design variable, one line per estimator
#'
#' @param data Cell summaries from summarise_sim().
#' @param y Unquoted column for the y axis.
#' @param ylab Y axis label.
#' @param x Unquoted design column for the x axis.
#' @param ests Estimators to include, in the canonical lab_est order.
#' @param from_zero Logical; extend the y axis to include 0.
#' @return A ggplot object; the caller adds scales and facets.
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
#' RMSE and |bias| side by side for one outcome
#'
#' @param data Cell summaries from summarise_sim().
#' @param oc Outcome name.
#' @param x Unquoted design column for the x axis.
#' @param ests Estimators to include.
#' @param extra A ggplot layer, or list of layers, appended to both panels.
#' @param title Plot title.
#' @return A patchwork object with one shared legend.
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

#' Share of MSE from squared bias and from variance, per estimator
#'
#' @param data Cell summaries from summarise_sim().
#' @param ests Estimators to include (rows of the facet grid).
#' @param x Unquoted design column for the x axis.
#' @param cols Unquoted column for the facet columns.
#' @param x_scale A ggplot x scale.
#' @param xlab X axis label.
#' @param col_lab A ggplot labeller for the facets.
#' @param title Plot title.
#' @return A ggplot object (stacked area plot).
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
#' Performance table: estimators by design columns
#'
#' Each cell shows the top metric with the bottom metric beneath it; rows are
#' ordered by the mean of order_by. With two design columns the first becomes
#' column spanners; with group, rows are grouped (for example by outcome).
#'
#' @param data Cell summaries from summarise_sim().
#' @param ests Estimators to include.
#' @param cols One or two design column names.
#' @param top Metric on top of each cell: "rmse", "bias", "abs_bias" or "empse".
#' @param bottom Metric beneath it, same choices.
#' @param digits Decimal places.
#' @param order_by Metric used to order the rows.
#' @param title Table title.
#' @param subtitle Subtitle; defaults to "<top> (<bottom>)".
#' @param group Optional column name for row groups.
#' @param span_lab Function turning a spanner value into its label.
#' @param span_desc Logical; order the spanners in decreasing order.
#' @return A gt table.
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
#' Share of replications whose path endpoint sits above lambda = .05
#'
#' @param data Long results or wide diagnostics.
#' @param by Grouping columns.
#' @param cols Which element of by becomes the table columns.
#' @param outcome_name Optional outcome to filter on.
#' @param title Optional title (markdown).
#' @return A gt table of percentages.
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

#' Table of misspecified minus correct RMSE and |bias|
#'
#' @param d Output of delta_ms().
#' @param ests Estimators to include.
#' @param at_n Sample size to show.
#' @param ov Overlap values to show.
#' @param digits Decimal places.
#' @param title Table title.
#' @return A gt table grouped by outcome, with overlap as the stub.
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
#' Delta RMSE and delta |bias| against overlap, faceted by n
#'
#' @param data Output of delta_ms().
#' @param oc Outcome name.
#' @param ests Estimators to include.
#' @param ov Overlap values to show.
#' @param title Plot title.
#' @return A patchwork object.
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

#' |Bias| and SD against overlap, correct versus misspecified
#'
#' @param data Cell summaries with a logical misspec column.
#' @param ests Estimators to include.
#' @param at_n Sample size to show.
#' @param ov Overlap values to show.
#' @return A ggplot object faceted by metric and by estimator x outcome.
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
#' Heatmap of extreme estimates (|tau_hat| > 10) for one estimator
#'
#' Tiles are coloured by the count per cell; labels show count / total over
#' the whole figure.
#'
#' @param data Cell summaries from summarise_sim() (uses n_extreme).
#' @param est Estimator name.
#' @param x Unquoted design column for the x axis.
#' @param xlab X axis label.
#' @param x_desc Logical; order the x levels in decreasing order.
#' @param ov Optional subset of x values.
#' @param title Plot title; defaults to the estimator label.
#' @return A ggplot object.
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
#' Median selected lambda, path endpoint and glmnet lambda.min against a design variable
#'
#' Optional q10-q90 ribbon; dashed black line at sqrt(log(p) / n) when rate is
#' TRUE (needs p in dg or p_dim); dotted lines at the fixed lambdas.
#'
#' @param dg Output of summarise_diag().
#' @param x Unquoted design column for the x axis (the rate line assumes n).
#' @param stats Which lam_ statistics to draw (see lam_meta).
#' @param ribbon Logical; draw the q10-q90 ribbon.
#' @param rate Logical; draw the sqrt(log(p) / n) reference.
#' @param ref_fixed Fixed lambdas to mark, or NULL.
#' @param p_dim Number of covariates when p is not a column of dg.
#' @return A ggplot object with a log10 y axis.
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

# attained max|SMD| vs its ceiling: med (q10-q90) by arm; hline at the fixed
# lambda for at = "05"/"10". For at = "cv"/"0" the bound varies per rep --
# read jointly with kkt_check / lambda_plot.
#' Attained max |SMD| against a design variable, by arm
#'
#' @param dg Output of summarise_diag().
#' @param x Unquoted design column for the x axis.
#' @param at Which weights: "05" or "10" (fixed lambda, drawn as a dotted
#'   ceiling), "cv" or "0" (path endpoint).
#' @param ribbon Logical; draw the q10-q90 ribbon.
#' @return A ggplot object.
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

#' Median active-set size against a design variable
#'
#' @param dg Output of summarise_diag() (uses the nnz_*_med columns).
#' @param x Unquoted design column for the x axis.
#' @return A ggplot object.
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

# WHERE: cell RMSE against the lambda each estimator sits at -- balnet0 at
# the median attained endpoint, .05/.10 fixed, balnetrate at sqrt(log p / n).
# balnetcv enters as a horizontal line, its selected-lambda q10-q90 mass as a
# band on the x axis: CV beats fixed lambda only if the line dips under the
# curve where its band sits. TODO once simulate_grid saves tau_path, replace
# these four points with the full tau(lambda) risk profile.
#' Cell RMSE against the lambda each fixed-lambda estimator sits at
#'
#' Points: balnet0 at the median attained endpoint, balnet05 and balnet10 at
#' their fixed lambdas, balnetrate at sqrt(log(p) / n). The CV estimator is a
#' horizontal line with its selected-lambda q10-q90 band on the x axis.
#'
#' @param data Long results from load_sim().
#' @param oc Outcome name.
#' @param facets Facet specification, as from vars().
#' @param labeller Facet labeller.
#' @return A ggplot object with a log10 x axis.
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
#' Mean CV estimate within quantile bins of the selected lambda
#'
#' Replications are binned on lam_cv within each cell; the binned mean of the
#' CV estimate is drawn with a +/- SD ribbon.
#'
#' @param data Long results from load_sim().
#' @param oc Outcome name.
#' @param bins Number of quantile bins.
#' @param ref_slope Optional slope of a dotted reference line through 0.
#' @param facets Facet specification, as from vars().
#' @param labeller Facet labeller.
#' @return A ggplot object.
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
#' Squared-error regret of the CV estimate against a fixed-lambda comparator
#'
#' Per replication, balnetcv^2 minus ref^2 (true tau = 0), binned on the
#' selected lambda within each cell; points are coloured by the sign of the
#' binned mean regret.
#'
#' @param data Long results from load_sim().
#' @param oc Outcome name.
#' @param ref Comparator estimator name.
#' @param bins Number of quantile bins.
#' @param facets Facet specification, as from vars().
#' @param labeller Facet labeller.
#' @return A ggplot object.
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

# spaghetti of per-rep CV loss curves; black points mark each rep's selected
# lambda. Filter to one arm / cell subset before plotting, facet as needed.
#' Spaghetti plot of per-replication CV loss curves
#'
#' @param cur Curves from load_cv_curves(), filtered to one arm and cell subset.
#' @param max_reps Number of replications to draw.
#' @param ref_fixed Lambdas marked with dotted vertical lines.
#' @return A ggplot object; black points mark each curve's minimum.
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

