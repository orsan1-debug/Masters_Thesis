# tables and figs from the verified summaries in output/thesis_summaries/
# nothing computed from results/ here (rbar.R does at build). load via handoff/load.R (own env)

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(gt)

#' x or y if x is NULL. base R only has it from 4.4
`%||%` <- function(x, y) if (is.null(x)) y else x

# labels and number formats

# this order = row order of every table
arm_labs <- c(balnet0 = "balnet", balnetcv = "balnet (CV)",
              abw_balnet0 = "ABW-balnet",
              glmnetcv_hajek = "glmnet (Hájek)", glmnetcv_ht = "glmnet (HT)",
              aipw_glmnetcv = "AIPW-glmnet", ols = "OLS",
              oracle_hajek = "Oracle (Hájek)", oracle_ht = "Oracle (HT)",
              oracle_aipw = "Oracle AIPW", unadjusted = "Unadjusted",
              abw_balnet0_ortrue = "ABW-balnet (true OR)",
              aipw_glmnetcv_ortrue = "AIPW-glmnet (true OR)")

outcome_models <- c("linear", "quad1", "exp")
outcome_labs <- c(linear = "Linear", quad1 = "Quadratic", exp = "Exponential")

# tables give the full spec name. figure panels are narrower
block_labs <- c(correct = "Correct specification",
                misspec = "Misspecified propensity score",
                both = "Both specifications")
block_labs_short <- c(correct = "Correct", misspec = "Misspecified")

font_serif <- c("Times New Roman", "Times", "serif")

#' plain decimals. no trailing zeros or sci notation
format_plain <- function(x) {
  format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
}

#' no leading zero (.05) like the share tiles. qmd passes it as x_fmt so keep the name
fmt_dot <- function(x) sub("^0\\.", ".", format_plain(x))

#' n in thousands for an axis titled n (thousands)
format_thousands <- function(x) as.character(x / 1000)

#' n as 500 1k ... 400k
format_k <- function(x) {
  ifelse(x >= 1000, paste0(x / 1000, "k"), format_plain(x))
}

#' fixed decimals w/ a true minus sign. a hyphen does not typeset as one
format_signed <- function(x, digits = 3) {
  out <- sprintf(paste0("%.", digits, "f"), x)
  # no sign if it rounds to 0. it would claim more than the digits shown
  out <- sub("^-(0(\\.0+)?)$", "\\1", out)
  sub("^-", "−", out)
}

#' label a pair A vs B. a shared (true OR) suffix is written once at the end
pair_name <- function(A, B) {
  lab_A <- unname(arm_labs[A])
  lab_B <- unname(arm_labs[B])
  tag <- " (true OR)"
  shared <- endsWith(lab_A, tag) & endsWith(lab_B, tag)
  ifelse(shared,
         paste0(sub(tag, "", lab_A, fixed = TRUE), " vs ",
                sub(tag, "", lab_B, fixed = TRUE), tag),
         paste(lab_A, "vs", lab_B))
}

#' display rank. pairs of one A together in the order A first shows up in pairs
pair_rank <- function(A, B, pairs) {
  position <- match(paste(A, B), paste(pairs$A, pairs$B))
  match(A, unique(pairs$A)) * (nrow(pairs) + 1) + position
}

#' booktabs look for a gt table. rules at top and under labels and at bottom only
style_booktabs <- function(g, font = font_serif) {
  g <- g |>
    opt_table_lines("none") |>
    tab_options(
      table.border.top.style = "solid", table.border.top.width = px(1.5),
      table.border.top.color = "black",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = px(0.75),
      column_labels.border.bottom.color = "black",
      table.border.bottom.style = "solid",
      table.border.bottom.width = px(1.5),
      table.border.bottom.color = "black",
      table.font.names = font,
      # pdf tables at 8pt for the wide appendix tables. html keeps gt default
      table.font.size = if (knitr::is_latex_output()) "8pt" else "16px"
    )
  if (!isTRUE(knitr::is_latex_output())) return(g)
  # a table over a page cannot be a float (pdf render loops) so only short ones set whole
  # long ones break across pages w/ column labels repeated
  long <- nrow(g[["_data"]]) + length(g[["_row_groups"]]) > latex_long_lines
  tab_options(g, latex.use_longtable = long, latex.header_repeat = long,
              latex.tbl.pos = "H")
}

# rows + row groups after which a pdf table breaks. about this many fit at 8pt
latex_long_lines <- 55

# reading the summaries

#' load one stored summary. warns if a check failed or a results file changed since build
read_summary <- function(id) {
  path <- here::here("output", "thesis_summaries", paste0(id, ".rds"))
  if (!file.exists(path)) {
    stop("no summary '", id, "': run handoff/scripts/build_summaries.R")
  }
  s <- readRDS(path)
  failed <- c(s$checks$check[!s$checks$ok],
              s$verify$check[s$verify$level == "error" & !s$verify$ok])
  if (is.null(s$verify)) failed <- c(failed, "summary stored without verify")
  if (length(failed)) {
    warning("summary '", id, "' failed: ", paste(failed, collapse = "; "))
  }
  # content decides not mtime. copying files changes mtime so md5 when size same
  # a missing results file is not a change. thesis must render without them
  built_from <- s$info$files
  results <- here::here("results", s$info$dir %||% "aipw", built_from$file)
  size <- file.size(results)
  changed <- !is.na(size) & size != built_from$size
  touched <- !is.na(size) & !changed &
    file.mtime(results) != built_from$mtime
  if (any(touched)) {
    changed[touched] <- unname(tools::md5sum(results[touched])) !=
      built_from$md5[touched]
  }
  if (any(changed)) {
    warning("results changed since summary '", id, "' was built (",
            paste(built_from$file[changed], collapse = ", "),
            "): rerun handoff/scripts/build_summaries.R")
  }
  s
}

# tables of averages
# A = estimator B = comparator. R_c = RMSE_A / RMSE_B and Rbar = 100 {1 - exp(mean log R_c)}

#' Rbar(Theta) by outcome and sbar_m(Theta) of A and B. one row per pair (atbl-q1-adv etc)
table_avg <- function(s, set = "all", blocks = NULL, pairs = NULL,
                      note = NULL, cell_counts = TRUE, title = NULL,
                      by_outcome = TRUE, shares = TRUE,
                      theta = c("pair", "common"), font = font_serif) {
  theta <- match.arg(theta)
  # theta pair = each pair over the cells it keeps. common = cells every estimator keeps
  blocks <- blocks %||% unique(s$avg$block)
  pairs <- pairs %||% s$info$pairs
  grouped <- length(blocks) > 1
  avg <- s$avg |>
    filter(.data$set == .env$set, block %in% blocks,
           .data$theta == .env$theta) |>
    semi_join(pairs, by = c("A", "B"))
  outcomes_shown <- intersect(outcome_models, avg$outcome)
  pooled <- filter(avg, outcome == "pooled")
  rbar_cols <- paste0("R_", c("pooled", if (by_outcome) outcomes_shown))
  # count goes last next to the avgs it qualifies and only if a row lost cells
  # under common every row loses the same cells so one note sentence says it
  show_left_out <- cell_counts && theta == "pair" &&
    any(pooled$n_dropped > 0)

  rows <- avg |>
    select(block, A, B, outcome, R_bar) |>
    pivot_wider(names_from = outcome, values_from = R_bar,
                names_prefix = "R_") |>
    left_join(select(pooled, block, A, B, s_A = s_bar_A, s_B = s_bar_B,
                     left_out = n_dropped),
              by = c("block", "A", "B")) |>
    arrange(match(block, blocks), pair_rank(A, B, pairs)) |>
    # estimator named on its first row only so its pairs read as one
    mutate(est_A = if_else(duplicated(paste(block, A)), "",
                           unname(arm_labs[A])),
           est_B = unname(arm_labs[B]),
           block = unname(block_labs[as.character(block)])) |>
    select(all_of(c("block", "est_A", "est_B", rbar_cols,
                    if (shares) c("s_A", "s_B"),
                    if (show_left_out) "left_out")))

  g <- if (grouped) {
    gt(rows, groupname_col = "block") |>
      tab_style(cell_text(weight = "bold"), cells_row_groups())
  } else {
    gt(select(rows, -block))
  }
  if (!is.null(title)) {
    g <- g |>
      tab_header(title = title) |>
      opt_align_table_header("left")
  }
  share_labels <- if (grouped) list(s_A = "A", s_B = "B") else
    list(s_A = md("$m = A$"), s_B = md("$m = B$"))
  if (shares) {
    g <- g |>
      tab_spanner(md("$\\bar s_m(\\Theta)$"), c(s_A, s_B)) |>
      cols_label(.list = share_labels) |>
      fmt_number(c(s_A, s_B), decimals = 2)
  }

  g |>
    tab_spanner(md("$\\bar R(\\Theta)$"), all_of(rbar_cols)) |>
    cols_label(.list = c(
      list(est_A = "Estimator (A)", est_B = "Comparator (B)",
           R_pooled = "All"),
      if (by_outcome) set_names(as.list(unname(outcome_labs[outcomes_shown])),
                                paste0("R_", outcomes_shown)),
      if (show_left_out) list(left_out = "Cells left out")
    )) |>
    fmt_number(all_of(rbar_cols), decimals = 1) |>
    cols_align("left", c(est_A, est_B)) |>
    cols_align("right", -c(est_A, est_B)) |>
    style_booktabs(font) |>
    # hairline above each estimator's first row to split the estimators
    tab_style(cell_borders(sides = "top", color = "#BDBDBD",
                           weight = px(0.5)),
              cells_body(rows = est_A != "" & seq_along(est_A) > 1)) |>
    tab_source_note(md(paste(c(
      if (cell_counts) note_cells(pooled, length(outcomes_shown)),
      note_avg_mcse(avg, shares),
      if (cell_counts) note_dropped(pooled, theta),
      note
    ), collapse = " ")))
}

#' note on cells averaged. only when every row has the same cells (no per row count)
note_cells <- function(pooled, n_outcomes) {
  if (n_distinct(pooled$n_cells) > 1) return(NULL)
  k <- pooled$n_cells[1]
  sprintf("%d cells, %d per outcome.", k, k / n_outcomes)
}

#' note w/ largest MCSE of Rbar (points) and sbar so the reader can judge the decimals
note_avg_mcse <- function(avg, shares = TRUE) {
  if (!shares) {
    return(sprintf("Largest MCSE: %.1f points ($\\bar R(\\Theta)$).",
                   max(avg$R_bar_mcse, na.rm = TRUE)))
  }
  sprintf(paste("Largest MCSE: %.1f points ($\\bar R(\\Theta)$),",
                "%.3f ($\\bar s_m(\\Theta)$)."),
          max(avg$R_bar_mcse, na.rm = TRUE),
          max(avg$s_bar_A_mcse, avg$s_bar_B_mcse, na.rm = TRUE))
}

#' note on cells lost to missing ests. pair defines the left out column. common = one count
note_dropped <- function(pooled, theta = "pair") {
  lost <- pooled$n_dropped[pooled$n_dropped > 0]
  if (!length(lost)) return(NULL)
  if (theta == "common") {
    return(sprintf(paste("Every pair over the same cells: %d left out where",
                         "an arm lacks more than 1%% of its estimates."),
                   max(lost)))
  }
  paste("Cells left out: where A or B lacks more than 1% of its estimates;",
        "the pair's averages skip them.")
}

#' Rbar per pair and sbar of one pair at each level of some by-sets
table_sets <- function(s, sets, pairs, share_pair = pairs[1, ],
                       factor_labs = NULL, block = s$info$specs[1],
                       note = NULL, theta = c("pair", "common"),
                       font = font_serif) {
  theta <- match.arg(theta)
  factor_labs <- factor_labs %||% character()
  avg <- s$avg |>
    filter(.data$block == .env$block, set %in% names(sets),
           outcome == "pooled", .data$theta == .env$theta) |>
    semi_join(pairs, by = c("A", "B")) |>
    mutate(pair = paste(A, B),
           factor_lab = coalesce(unname(factor_labs[factor]), factor),
           group = case_when(
             is.na(factor) ~ unname(sets[set]),
             unname(sets[set]) == "" ~ sub("^(.)", "\\U\\1", factor_lab,
                                           perl = TRUE),
             .default = paste0(sets[set], ", ", factor_lab)),
           value = coalesce(value, "all cells"),
           set = factor(set, names(sets)))
  shares <- avg |>
    filter(A == share_pair$A, B == share_pair$B) |>
    select(set, factor, value, s_A = s_bar_A, s_B = s_bar_B)
  rbar_cols <- paste(pairs$A, pairs$B)
  rows <- avg |>
    select(set, group, factor, value, pair, R_bar) |>
    pivot_wider(names_from = pair, values_from = R_bar) |>
    left_join(shares, by = c("set", "factor", "value")) |>
    arrange(set) |>
    select(group, value, all_of(rbar_cols), s_A, s_B)

  gt(rows, groupname_col = "group", rowname_col = "value") |>
    cols_label(.list = c(
      set_names(as.list(pair_name(pairs$A, pairs$B)), rbar_cols),
      list(s_A = unname(arm_labs[share_pair$A]),
           s_B = unname(arm_labs[share_pair$B]))
    )) |>
    tab_spanner(md("$\\bar R(\\Theta)$, A vs B"), all_of(rbar_cols)) |>
    tab_spanner(md("$\\bar s_m(\\Theta)$"), c(s_A, s_B)) |>
    fmt_number(all_of(rbar_cols), decimals = 1) |>
    fmt_number(c(s_A, s_B), decimals = 2) |>
    style_booktabs(font) |>
    tab_source_note(md(paste(c(note_avg_mcse(avg), note), collapse = " ")))
}

# tables per cell

#' RMSE (bias) per estimator and cell (appendix). dash = over 1% ests missing. blank = not run
table_cells <- function(s, block = s$info$specs[1], x = s$info$x,
                        by = s$info$by, arms = NULL,
                        layout = c("wide", "long"), row_cols = NULL,
                        design_labs = NULL, x_fmt = NULL, x_name = x,
                        x_title = NULL, by_fmt = NULL, title = NULL,
                        note = NULL, digits = 3, whole_from = Inf,
                        font = font_serif) {
  layout <- match.arg(layout)
  x_fmt <- x_fmt %||% if (x == "n") format_thousands else format_plain
  x_title <- x_title %||% if (x == "n") "*n* (thousands)" else x_name
  by_fmt <- by_fmt %||% \(v) paste(by, "=", v)
  design_label <- \(v) coalesce(unname((design_labs %||% character())[v]), v)
  perf <- filter(s$perf, spec == block)
  arms <- intersect(names(arm_labs), arms %||% unique(perf$estimator))
  perf <- filter(perf, estimator %in% arms)
  # values >= whole_from print to the unit. MCSE too big for decimals (glmnet HT tails)
  shown <- \(v) if_else(abs(v) >= whole_from, format_signed(v, 0),
                        format_signed(v, digits))
  arm_cells <- perf |>
    mutate(cell = if_else(dropped, "—",
                          paste0(shown(rmse), " (", shown(bias), ")")),
           Outcome = factor(outcome_labs[outcome],
                            outcome_labs[intersect(outcome_models, outcome)]))

  if (layout == "wide") {
    levels_x <- arm_cells |>
      distinct(across(all_of(c(by, x)))) |>
      arrange(across(all_of(c(by, x))))
    # col ids carry the level of by so the spanners can pick them
    column_id <- \(d) if (is.null(by)) paste0("x_", d[[x]]) else
      paste0("b_", d[[by]], "|x_", d[[x]])
    cols <- column_id(levels_x)
    rows <- arm_cells |>
      mutate(col = column_id(arm_cells),
             Arm = factor(unname(arm_labs[estimator]),
                          unname(arm_labs[arms]))) |>
      select(Outcome, Arm, col, cell) |>
      pivot_wider(names_from = col, values_from = cell) |>
      arrange(Outcome, Arm) |>
      select(Outcome, Arm, all_of(cols))
    g <- gt(rows, groupname_col = "Outcome", rowname_col = "Arm") |>
      cols_label(.list = set_names(as.list(x_fmt(levels_x[[x]])), cols)) |>
      tab_stubhead(md(x_title))
    if (!is.null(by)) {
      for (v in unique(levels_x[[by]])) {
        g <- tab_spanner(g, by_fmt(v), starts_with(paste0("b_", v, "|")))
      }
    }
    mcse_cols <- c(by, x)
    mcse_col_names <- c(by, x_name)
  } else {
    row_cols <- row_cols %||% c(by, x)
    # sort on the numbers before printing so 10 sorts after 2
    rows <- arm_cells |>
      arrange(Outcome, across(all_of(row_cols))) |>
      mutate(Arm = unname(arm_labs[estimator]),
             across(all_of(row_cols), format_plain)) |>
      select(Outcome, all_of(row_cols), Arm, cell) |>
      pivot_wider(names_from = Arm, values_from = cell) |>
      select(Outcome, all_of(row_cols), any_of(unname(arm_labs[arms])))
    g <- gt(rows, groupname_col = "Outcome") |>
      cols_label(.list = set_names(as.list(design_label(row_cols)),
                                   row_cols)) |>
      cols_align("left", all_of(row_cols))
    mcse_cols <- row_cols
    mcse_col_names <- design_label(row_cols)
  }
  if (!is.null(title)) {
    g <- g |>
      tab_header(title = title) |>
      opt_align_table_header("left")
  }
  g |>
    cols_align("right", any_of(unname(arm_labs[arms]))) |>
    cols_align("right", starts_with(c("x_", "b_"))) |>
    style_booktabs(font) |>
    tab_style(cell_text(align = "left", whitespace = "nowrap"),
              cells_stub()) |>
    sub_missing(everything(), missing_text = "") |>
    tab_source_note(md(paste(c(
      note_cell_mcse(perf, mcse_cols, mcse_col_names),
      note_missing(perf, s$info$miss_max),
      if (any(abs(c(perf$rmse, perf$bias)) >= whole_from, na.rm = TRUE)) {
        sprintf("Values of %s or more are shown to the unit.",
                format_plain(whole_from))
      },
      note
    ), collapse = " ")))
}

#' pdf only. fixed col widths (pt incl padding) so long headers wrap and the table fits
pdf_col_width <- function(g, cols, width) {
  if (!isTRUE(knitr::is_latex_output())) return(g)
  # gt evaluates the formula away from here so col names and width go in as values
  w <- rlang::new_formula(rlang::expr(c(!!!rlang::syms(cols))),
                          rlang::expr(px(!!(width / 0.75))))
  cols_width(g, .list = list(w))
}

#' note w/ largest MCSE of RMSE and of bias and the cell each comes from
note_cell_mcse <- function(perf, cols, col_names = cols) {
  kept <- filter(perf, !dropped)
  located <- \(row) paste(sprintf("%s = %s", col_names,
                                  map_chr(cols, \(v) format_plain(row[[v]]))),
                          collapse = ", ")
  r <- kept[which.max(kept$mcse_rmse), ]
  b <- kept[which.max(kept$mcse_bias), ]
  sprintf(
    "Largest MCSE: %s (RMSE; %s, %s, %s) and %s (bias; %s, %s, %s).",
    format_plain(signif(r$mcse_rmse, 2)), arm_labs[r$estimator],
    tolower(outcome_labs[r$outcome]), located(r),
    format_plain(signif(b$mcse_bias, 2)), arm_labs[b$estimator],
    tolower(outcome_labs[b$outcome]), located(b)
  )
}

#' note on missing ests per estimator. reps and cells w/ any and cells over miss_max (dash)
note_missing <- function(perf, miss_max = 0.01) {
  # same reps missing in every outcome model (true of every stored cell) so count one
  # dplyr:: since data.table masks first() when attached
  m <- perf |>
    filter(missing > 0) |>
    summarise(cells = n_distinct(level), over = n_distinct(level[dropped]),
              reps = sum(missing[outcome == dplyr::first(outcome)]),
              .by = estimator) |>
    arrange(match(estimator, names(arm_labs)))
  if (!nrow(m)) return(NULL)
  paste0("Missing estimates: ",
         paste(sprintf("%s, %d in %d of %d cells (%d above %s%%, shown as —)",
                       arm_labs[m$estimator], m$reps, m$cells,
                       n_distinct(perf$level), m$over,
                       format_plain(100 * miss_max)),
               collapse = "; "),
         ".")
}

#' IPW-balnet linear bias as % of the imbalance bound once the path truncates
table_sp_bound <- function(s, font = font_serif) {
  bound <- s$extra$bound
  bound |>
    transmute(s, decay, corr, block = paste0("n", n, "_pn", pn), pct) |>
    pivot_wider(names_from = block, values_from = pct) |>
    arrange(s, desc(decay), corr) |>
    gt() |>
    cols_label(s = md("*s*"), decay = md("*ν*"), corr = "Covariates",
               n500_pn0.5 = "0.5", n500_pn1 = "1",
               n1000_pn0.5 = "0.5", n1000_pn1 = "1") |>
    tab_spanner(md("*n* = 500, *p*/*n*"), c(n500_pn0.5, n500_pn1)) |>
    tab_spanner(md("*n* = 1000, *p*/*n*"), c(n1000_pn0.5, n1000_pn1)) |>
    fmt_number(c(n500_pn0.5, n500_pn1, n1000_pn0.5, n1000_pn1),
               decimals = 0) |>
    style_booktabs(font) |>
    tab_source_note(md(sprintf(paste(
      "Largest MCSE: %.1f points; values above 100 are within 1 MCSE of 100.",
      "The $\\ell_1$ norm is $\\sum_{j \\le s} j^{-\\nu} / \\sigma_g$."),
      max(bound$pct_mcse))))
}

# heat maps

# plotmath symbols for panel titles and legends in thesis letters
symbol_ratio <- quote(100 * (1 - italic(R)[italic(c)]))
symbol_share <- quote(italic(s)[italic(m) * "," * italic(c)])

#' heat map theme. colour bar shrinks w/ the rows so a short map is not stretched
theme_heat <- function(n_rows, family = "serif", base_size = 9,
                       bar_per_row = 0.22) {
  bar_height <- min(1, max(0.6, bar_per_row * n_rows))
  theme_minimal(base_size = base_size, base_family = family) +
    theme(panel.grid = element_blank(),
          panel.spacing.x = unit(0.8, "lines"),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(size = base_size - 2),
          legend.position = "right",
          legend.key.width = unit(0.2, "cm"),
          legend.key.height = unit(bar_height / 5, "in"),  # bar = 5 keys
          legend.title = element_text(size = base_size - 1),
          legend.text = element_text(size = base_size - 2),
          legend.margin = margin(0, 0, 0, 2),
          plot.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = base_size),
          strip.text.y = element_text(angle = 0, hjust = 0),
          plot.title = element_text(size = base_size + 1),
          plot.subtitle = element_text(size = base_size, face = "italic"),
          plot.title.position = "plot")
}

#' colour limits on log2 R_c. one number is symmetric (2 = x4 either way) or give lo and hi
heat_fill_limits <- function(lim) {
  lims <- if (length(lim) == 1) c(-lim, lim) else sort(lim)
  stopifnot(length(lims) == 2, lims[1] <= 0, lims[2] >= 0)
  lims
}

#' legend breaks as ratios (powers of 2 inside the limits). coarse = 0.5 1 2 for short bars
heat_legend_breaks <- function(lims, coarse = FALSE) {
  k <- if (coarse) -1:1 else seq(floor(lims[1]), ceiling(lims[2]))
  2^k[k >= lims[1] & k <= lims[2]]
}

#' heat map facets. outcomes across and a row of panels per level of by
heat_facets <- function(by, by_fmt) {
  if (is.null(by)) return(facet_wrap(~outcome))
  facet_grid(rows = vars(.data[[by]]), cols = vars(outcome),
             labeller = labeller(.rows = as_labeller(by_fmt),
                                 .cols = label_value))
}

#' heat map title. the axes and strips name what varies
heat_title <- function(what, symbol) {
  bquote(.(what) * " per cell, " * .(symbol))
}

#' order the outcome models and x levels for the tiles
heat_tile_axes <- function(d, x, x_levels, x_fmt) {
  d |>
    mutate(outcome = factor(outcome_labs[outcome],
                            outcome_labs[intersect(outcome_models, outcome)]),
           x_f = factor(.data[[x]], x_levels, x_fmt(x_levels)))
}

#' panel (a) 100(1 - R_c) per pair and cell. log2 fill so x2 either way is equally far from white
heat_panel_ratio <- function(cells, pair_levels, by, by_fmt, title, x_lab,
                             spec, flag, lims, family, base_size, text_size) {
  n_rows <- length(pair_levels)
  ratio_breaks <- heat_legend_breaks(lims, coarse = n_rows < 3)
  tiles <- cells |>
    mutate(pair = factor(pair_name(A, B), rev(pair_levels)),
           pct = round(100 * (1 - R_c)) + 0,   # + 0 so no tile reads -0
           label = if_else(is.na(pct), "—",
                           sub("-", "−", sprintf("%.0f", pct))),
           face = if (is.null(flag)) "plain" else
             if_else(!is.na(pct) & abs(log(R_c)) < flag * log_R_c_mcse,
                     "italic", "plain"),
           # white text where the fill is dark
           ink = if_else(!is.na(pct) &
                           abs(log2(R_c)) > 0.55 * max(abs(lims)),
                         "white", "grey15"))
  ggplot(tiles, aes(x_f, pair, fill = log2(R_c))) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = label, colour = ink, fontface = face),
              size = text_size, family = family) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F0EFEC", high = "#B2182B", midpoint = 0,
      limits = lims, oob = scales::squish, na.value = "grey85",
      breaks = log2(ratio_breaks),
      labels = sub("-", "−", round(100 * (1 - ratio_breaks))),
      name = as.expression(symbol_ratio),
      # reversed so blue (A has lower RMSE) is on top
      guide = guide_colourbar(title.position = "top", title.hjust = 0,
                              reverse = TRUE)
    ) +
    heat_facets(by, by_fmt) +
    labs(x = x_lab, y = NULL, title = title, subtitle = spec) +
    theme_heat(n_rows, family, base_size)
}

#' panel (b) share of bias per estimator and cell on 0 to 1
heat_panel_share <- function(shares, arm_levels, n_rows, by, by_fmt, title,
                             x_lab, spec, family, base_size, text_size) {
  tiles <- shares |>
    mutate(arm = factor(arm_labs[arm], rev(arm_labs[arm_levels])),
           label = if_else(is.na(share), "—",
                           sub("^0", "", sprintf("%.2f", share))),
           ink = if_else(!is.na(share) & share > 0.55, "white", "grey15"))
  ggplot(tiles, aes(x_f, arm, fill = share)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = label, colour = ink), size = text_size,
              family = family) +
    scale_colour_identity() +
    scale_fill_gradient(low = "#F7F7F7", high = "#3F2A7A",
                        limits = c(0, 1), na.value = "grey85",
                        breaks = c(0, 0.5, 1), labels = c("0", ".5", "1"),
                        name = as.expression(symbol_share),
                        guide = guide_colourbar(title.position = "top",
                                                title.hjust = 0)) +
    heat_facets(by, by_fmt) +
    labs(x = x_lab, y = NULL, title = title, subtitle = spec) +
    # tiles half as tall as panel (a) (heat_share_height())
    theme_heat(n_rows, family, base_size, bar_per_row = 0.11)
}

#' per cell heat maps. (a) 100(1 - R_c) per pair (b) share of bias per estimator
#' 2 specs = a row of panels each on one colour scale so they compare tile by tile
heat_cells <- function(s, x = s$info$x, blocks = NULL, pairs = NULL,
                       keep = NULL, by = NULL, arms = NULL, x_lab = NULL,
                       x_fmt = NULL, by_fmt = NULL, spec = NULL, flag = 2,
                       lim = 2, family = "serif", base_size = 9,
                       text_size = 2.3) {
  blocks <- blocks %||% s$info$specs
  pairs <- pairs %||% s$info$pairs
  is_n <- x == "n"
  x_lab <- x_lab %||% if (is_n) expression(italic(n) * " (thousands)") else x
  x_fmt <- x_fmt %||% if (is_n) format_thousands else format_plain

  cells <- s$cells |>
    filter(spec %in% blocks) |>
    semi_join(pairs, by = c("A", "B")) |>
    arrange(match(spec, blocks),
            match(paste(A, B), paste(pairs$A, pairs$B)))
  if (!is.null(keep)) cells <- cells[keep(cells), ]
  arm_levels <- intersect(names(arm_labs), c(cells$A, cells$B))
  if (!is.null(arms)) arm_levels <- intersect(arm_levels, arms)
  # panel (b) covers the cells of panel (a). a cell left out there is grey
  shares <- s$perf |>
    filter(estimator %in% arm_levels) |>
    semi_join(cells, by = c("spec", "level")) |>
    mutate(arm = estimator, share = if_else(dropped, NA_real_, share),
           mcse = if_else(dropped, NA_real_, mcse_share))

  if (length(blocks) > 1) {
    # spec goes in the facet rows on its own line so the strip stays narrow
    panel_of <- \(d) if (is.null(by)) block_labs_short[d$spec] else
      paste0(block_labs_short[d$spec], ",\n", by, " = ", d[[by]])
    panels <- if (is.null(by)) block_labs_short[blocks] else
      as.vector(t(outer(block_labs_short[blocks], sort(unique(cells[[by]])),
                        \(b, v) paste0(b, ",\n", by, " = ", v))))
    cells <- mutate(cells, panel = factor(panel_of(cells), panels))
    shares <- mutate(shares, panel = factor(panel_of(shares), panels))
    by <- "panel"
    by_fmt <- identity
  }
  by_fmt <- by_fmt %||% \(v) paste(by, "=", v)

  x_levels <- sort(unique(cells[[x]]))
  cells <- heat_tile_axes(cells, x, x_levels, x_fmt)
  shares <- heat_tile_axes(shares, x, x_levels, x_fmt)
  pair_levels <- unique(pair_name(cells$A, cells$B)[
    order(pair_rank(cells$A, cells$B, pairs))])
  n_panels <- if (is.null(by)) 1 else length(unique(cells[[by]]))

  list(
    rmse = heat_panel_ratio(cells, pair_levels, by, by_fmt,
                            heat_title("Relative RMSE", symbol_ratio),
                            x_lab, spec, flag, heat_fill_limits(lim), family,
                            base_size, text_size),
    share = heat_panel_share(shares, arm_levels,
                             length(arm_levels) * n_panels, by, by_fmt,
                             heat_title("Share of bias", symbol_share),
                             x_lab, spec, family, base_size, text_size),
    height = heat_height(cells, by),
    note = heat_note(cells, flag, lim),
    share_height = heat_share_height(length(arm_levels), n_panels),
    share_note = heat_note_share(shares$mcse, shares$arm,
                                 anyNA(shares$share))
  )
}

#' height of a share heat map drawn alone. rows less tall than wide so it reads as strips
heat_share_height <- function(rows, panels = 1, per_row = 0.13, extra = 1) {
  rows * panels * per_row + extra + 0.25 * (panels - 1)
}

#' caption for a share heat map. largest MCSE (approx if mcse_approx_arm) and grey cells
heat_note_share <- function(mcse, arm, grey = FALSE) {
  m_max <- max_where(mcse)
  m_rest <- max_where(mcse, arm != mcse_approx_arm)
  txt <- sprintf("The largest MCSE is %.3f.", m_max)
  if (m_max > m_rest && is.finite(m_rest)) {
    txt <- paste(txt, sprintf("It comes from %s, whose MCSE is approximate.",
                              arm_labs[[mcse_approx_arm]]))
  }
  paste(c(txt, if (grey) "Grey cells are left out for missing estimates."),
        collapse = " ")
}

#' fig height that keeps the tiles the same size whatever the number of rows
heat_height <- function(cells, by = NULL, per_row = 0.26, extra = 1) {
  rows <- max(nrow(distinct(cells, A, B)),
              length(unique(c(cells$A, cells$B))))
  panels <- if (is.null(by)) 1 else length(unique(cells[[by]]))
  rows * panels * per_row + extra + 0.25 * (panels - 1)
}

#' caption for heat_cells(). italic rule and largest MCSE per panel and range if one sided
heat_note <- function(cells, flag = 2, lim = 2) {
  lims <- heat_fill_limits(lim)
  pct <- \(k) sub("-", "−", round(100 * (1 - 2^k)))
  # symmetric scale reads off the legend. a one sided one needs its range
  scale_txt <- if (lims[1] != -lims[2]) {
    sprintf("The colour in (a) runs from %s to %s per cent lower RMSE.",
            pct(lims[2]), pct(lims[1]))
  }
  paste(c(
    if (!is.null(flag)) {
      sprintf("Italic values are within %g MCSE of no difference.", flag)
    },
    heat_note_mcse(cells),
    scale_txt,
    if (anyNA(cells$R_c)) "Grey cells are left out for missing estimates."
  ), collapse = " ")
}

# captions call this arm's MCSE approx. heavy tails (1000s of errors over 10 under misspec)
# put the delta method MCSE well below the bootstrap sd
mcse_approx_arm <- "glmnetcv_ht"

#' max of v over rows sel or -Inf if none so a comparison still works
max_where <- function(v, sel = TRUE) {
  keep <- sel & !is.na(v)
  if (any(keep)) max(v[keep]) else -Inf
}

#' MCSE sentences of heat_note(). says so when mcse_approx_arm sets the max
heat_note_mcse <- function(cells) {
  a <- 100 * cells$R_c * cells$log_R_c_mcse
  a_tail <- cells$A == mcse_approx_arm | cells$B == mcse_approx_arm
  b <- c(cells$s_A_mcse, cells$s_B_mcse)
  b_tail <- c(cells$A, cells$B) == mcse_approx_arm
  a_max <- max_where(a)
  b_max <- max_where(b)
  a_rest <- max_where(a, !a_tail)
  b_rest <- max_where(b, !b_tail)
  # flag only if the arm sets the max and some other estimator is shown
  a_flag <- a_max > a_rest && is.finite(a_rest)
  b_flag <- b_max > b_rest && is.finite(b_rest)
  largest <- sprintf("The largest MCSE is %.1f points in (a) and %.3f in (b).",
                     a_max, b_max)
  from <- sprintf("from %s, whose MCSE is approximate.",
                  arm_labs[[mcse_approx_arm]])
  if (a_flag && b_flag) {
    paste(largest, "Both come", from)
  } else if (a_flag) {
    paste(largest, "The first comes", from)
  } else if (b_flag) {
    paste(largest, "The second comes", from)
  } else {
    largest
  }
}

#' heat map of sbar_m(Theta) per estimator at each level of one factor (fig-ov-adv)
heat_share_avg <- function(s, set, x, blocks = c("correct", "misspec"),
                           pairs = NULL, arms = NULL, x_lab = x,
                           x_fmt = format_plain,
                           theta = c("pair", "common"),
                           family = "serif", base_size = 9,
                           text_size = 2.3) {
  theta <- match.arg(theta)
  pairs <- pairs %||% s$info$pairs
  avg <- s$avg |>
    filter(.data$set == .env$set, block %in% blocks, outcome != "pooled",
           .data$theta == .env$theta) |>
    semi_join(pairs, by = c("A", "B")) |>
    mutate(across(c(s_bar_A, s_bar_B, s_bar_A_mcse, s_bar_B_mcse),
                  \(v) if_else(n_dropped > 0, NA_real_, v)),
           !!x := as.numeric(value))
  arm_levels <- intersect(names(arm_labs), c(avg$A, avg$B))
  if (!is.null(arms)) arm_levels <- intersect(arm_levels, arms)
  # share is the same in each pair that kept every cell so take the first
  shares <- bind_rows(
    transmute(avg, block, outcome, across(all_of(x)), arm = A,
              share = s_bar_A),
    transmute(avg, block, outcome, across(all_of(x)), arm = B,
              share = s_bar_B)
  ) |>
    summarise(share = dplyr::first(na.omit(c(share, NA_real_))),
              .by = c(block, outcome, all_of(x), arm)) |>
    filter(arm %in% arm_levels) |>
    mutate(block = factor(block, blocks)) |>
    heat_tile_axes(x, sort(unique(avg[[x]])), x_fmt)

  p <- heat_panel_share(shares, arm_levels,
                        length(arm_levels) * length(blocks), "block",
                        \(v) unname(block_labs_short[as.character(v)]),
                        title = expression("Share of bias averaged over " *
                                             italic(n) * ", " *
                                             bar(italic(s))[italic(m)]),
                        x_lab = x_lab, spec = NULL, family = family,
                        base_size = base_size, text_size = text_size) +
    guides(fill = guide_colourbar(title = expression(bar(italic(s))[
      italic(m)]), title.position = "top", title.hjust = 0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  m_arm <- c(avg$A, avg$B)
  m <- c(avg$s_bar_A_mcse, avg$s_bar_B_mcse)[m_arm %in% arm_levels]
  list(
    share = p,
    # angled labels need more room
    height = heat_share_height(length(arm_levels), length(blocks),
                               extra = 1.3),
    note = heat_note_share(m, m_arm[m_arm %in% arm_levels],
                           anyNA(shares$share))
  )
}

# n x c grid. RMSE vs overlap (afig-app-nc)

# legend order. same marks as the strip figs so an estimator looks the same in every line fig
nc_arms <- c("balnet0", "abw_balnet0", "glmnetcv_hajek", "aipw_glmnetcv")

#' bands where at least half the balnet paths truncate. edges halfway between levels on log scale
nc_truncation_bands <- function(trunc) {
  trunc |>
    arrange(spec, n, overlap) |>
    group_by(spec, n) |>
    mutate(xmin = sqrt(overlap * lag(overlap, default = 0)),
           xmax = sqrt(overlap * lead(overlap, default = Inf)),
           run = cumsum(trunc >= 0.5 & lag(trunc < 0.5, default = TRUE))) |>
    filter(trunc >= 0.5) |>
    group_by(spec, n, run) |>
    summarise(xmin = min(xmin), xmax = max(xmax), .groups = "drop")
}

#' RMSE vs overlap c at some n of the n x c grid. both specs. cells w/ over 1% missing not drawn
#' bands +/- 1.96 MCSE (Morris et al. 2019 sec. 6.2). shaded where half the paths truncate
fig_nc <- function(s, n_sl, shade_label = "balnet truncates", title = NULL) {
  perf <- filter(s$perf, n %in% n_sl, estimator %in% nc_arms, !dropped)
  absent <- setdiff(n_sl, perf$n)
  if (length(absent)) stop("no cells at n = ", paste(absent, collapse = ", "))
  c_levels <- sort(unique(s$perf$overlap))
  # plotmath strips so n is italic as in the text
  spec_lab <- \(sp) factor(if_else(sp == "misspec", "Misspecified", "Correct"),
                           c("Correct", "Misspecified"))
  n_lab <- \(n) factor(paste0("italic(n) == ", n),
                       paste0("italic(n) == ", n_sl))
  d <- perf |>
    mutate(lo = pmax(rmse - 1.96 * mcse_rmse, 0),
           hi = rmse + 1.96 * mcse_rmse,
           spec = spec_lab(spec), n_lab = n_lab(n),
           outcome = factor(outcome_labs[outcome], outcome_labs),
           estimator = factor(estimator, nc_arms))
  bands <- nc_truncation_bands(filter(s$extra$trunc, n %in% n_sl)) |>
    mutate(spec = spec_lab(spec), n_lab = n_lab(n))
  corner <- filter(bands, spec == "Misspecified", n == max(n_sl)) |>
    slice_min(xmin, n = 1, with_ties = FALSE) |>
    mutate(outcome = factor(levels(d$outcome)[1], levels(d$outcome)),
           x = c_levels[1], label = shade_label)
  if (is.null(shade_label)) corner <- corner[0, ]
  ggplot(d, aes(overlap, rmse, colour = estimator, shape = estimator,
                linetype = estimator, group = estimator)) +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax,
                                ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "#EDE7DB") +
    geom_text(data = corner, aes(x = x, y = -Inf, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = -0.7, size = 2.6,
              family = "serif", colour = "grey30") +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = estimator), colour = NA,
                alpha = 0.18, show.legend = FALSE) +
    geom_line(linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 0.9) +
    ggh4x::facet_nested(outcome ~ spec + n_lab, scales = "free_y",
                        labeller = label_parsed) +
    # tick at every grid level. labels at the decades only
    scale_x_log10(breaks = c_levels,
                  labels = if_else(c_levels %in% c(0.01, 0.1, 1),
                                   fmt_dot(c_levels), "")) +
    scale_y_continuous(labels = fmt_ticks) +
    scale_colour_manual(values = npath_strip_pal[nc_arms],
                        labels = arm_labs[nc_arms]) +
    scale_fill_manual(values = npath_strip_pal[nc_arms], guide = "none") +
    scale_linetype_manual(values = npath_strip_lty[nc_arms], guide = "none") +
    scale_shape_manual(values = npath_strip_shape[nc_arms],
                       labels = arm_labs[nc_arms]) +
    guides(colour = guide_legend(override.aes = list(size = 1.8))) +
    expand_limits(y = 0) +
    labs(x = expression(paste("Overlap ", italic(c))), y = "RMSE",
         colour = NULL, shape = NULL, title = title) +
    theme_bw(base_size = 9, base_family = "serif") +
    theme(plot.title = element_text(size = 10, margin = margin(b = 3)),
          plot.title.position = "plot",
          legend.position = "bottom",
          legend.box.spacing = unit(0, "pt"),
          legend.key.spacing.x = unit(6, "pt"),
          legend.text = element_text(size = 7.5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2),
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          panel.spacing.x = unit(8, "pt"),
          axis.text = element_text(colour = "black"))
}

#' largest RMSE MCSE per estimator in fig_nc() for the caption. glmnet (Hájek) sets the max
nc_mcse <- function(s, n_sl) {
  drawn <- filter(s$perf, n %in% n_sl, estimator %in% nc_arms, !dropped)
  vapply(nc_arms, \(a) max(drawn$mcse_rmse[drawn$estimator == a]), 0)
}

# paths along one design factor. RMSE and |bias| (fig-n-misspec)

#' n as 250 1k 100k 1M for a log n axis
format_kM <- function(x) {
  ifelse(x >= 1e6, paste0(x / 1e6, "M"), format_k(x))
}

#' RMSE and |bias| rows w/ 95% MC bands (+/- 1.96 MCSE). |bias| band folded at 0
npath_rows <- function(perf) {
  rmse_rows <- mutate(perf, measure = "RMSE", y = rmse,
                      lo = pmax(rmse - 1.96 * mcse_rmse, 0),
                      hi = rmse + 1.96 * mcse_rmse)
  bias_rows <- mutate(perf, measure = "|Bias|", y = abs(bias),
                      b_lo = bias - 1.96 * mcse_bias,
                      b_hi = bias + 1.96 * mcse_bias,
                      lo = if_else(b_lo > 0 | b_hi < 0,
                                   pmin(abs(b_lo), abs(b_hi)), 0),
                      hi = pmax(abs(b_lo), abs(b_hi))) |>
    select(-b_lo, -b_hi)
  bind_rows(rmse_rows, bias_rows)
}

#' largest MCSE of RMSE and bias of the n path figs for the caption. faded arm and the rest
npath_mcse <- function(s, arms, block = s$info$specs[1], faded_arm = NULL) {
  drawn <- filter(s$perf, spec == block, estimator %in% arms, !dropped)
  faded <- drawn$estimator %in% faded_arm
  big <- \(v, sel) if (any(sel)) max(v[sel]) else NA_real_
  c(rmse_faded = big(drawn$mcse_rmse, faded),
    bias_faded = big(drawn$mcse_bias, faded),
    rmse_other = big(drawn$mcse_rmse, !faded),
    bias_other = big(drawn$mcse_bias, !faded))
}

#' max |diff| in RMSE and bias of 2 arms for x >= from. for a caption on why one is not drawn
npath_gap <- function(s, A, B, from, x = s$info$x, block = s$info$specs[1]) {
  d <- s$perf |>
    filter(spec == block, estimator %in% c(A, B), !dropped,
           .data[[x]] >= from) |>
    select(all_of(c(x, "outcome", "estimator", "rmse", "bias"))) |>
    pivot_wider(names_from = estimator, values_from = c(rmse, bias))
  gap <- \(m) max(abs(d[[paste0(m, "_", A)]] - d[[paste0(m, "_", B)]]),
                  na.rm = TRUE)
  c(rmse = gap("rmse"), bias = gap("bias"))
}

# paths along n in strips. RMSE and bias and share of bias (fig-n-misspec fig-n-refit)

# greys and one blue so figs read in print. line type and shape split shared inks
# ABW dark blue and balnet light. the 2 ABW arms never share a fig
npath_strip_pal <- c(balnet0 = "#56B4E9", abw_balnet0 = "#2166AC",
                     abw_balnet0_ortrue = "#2166AC",
                     glmnetcv_hajek = "#808080", aipw_glmnetcv = "#333333",
                     aipw_glmnetcv_ortrue = "#333333")
npath_strip_lty <- c(balnet0 = "solid", abw_balnet0 = "62",
                     abw_balnet0_ortrue = "62",
                     glmnetcv_hajek = "solid", aipw_glmnetcv = "22",
                     aipw_glmnetcv_ortrue = "42")
npath_strip_shape <- c(balnet0 = 16, abw_balnet0 = 1, abw_balnet0_ortrue = 1,
                       glmnetcv_hajek = 17, aipw_glmnetcv = 15,
                       aipw_glmnetcv_ortrue = 0)
npath_measures <- c("RMSE", "Bias", "Share of bias")

#' tick labels w/ one common number of decimals (min 1) so 0 1 2 never sits next to 0.0 0.5 1.0
fmt_ticks <- function(x) {
  ok <- !is.na(x)
  digits <- 1
  while (digits < 6 && any(abs(round(x[ok], digits) - x[ok]) > 1e-9)) {
    digits <- digits + 1
  }
  out <- rep(NA_character_, length(x))
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}

#' rows a strip fig draws. RMSE and |bias| as npath_rows() plus share w/ band kept in 0 to 1
npath_strip_rows <- function(s, arms, outcomes, measures, x, block,
                             separate) {
  perf <- filter(s$perf, spec == block, estimator %in% arms,
                 outcome %in% outcomes, !dropped)
  if (!nrow(perf)) stop("no cells of ", paste(arms, collapse = ", "))
  share <- mutate(perf, measure = "Share of bias", y = share,
                  lo = pmax(share - 1.96 * mcse_share, 0),
                  hi = pmin(share + 1.96 * mcse_share, 1))
  # signed so a path crossing zero reads as one drift not bias lost and regained
  signed <- mutate(perf, measure = "Bias", y = bias,
                   lo = bias - 1.96 * mcse_bias, hi = bias + 1.96 * mcse_bias)
  bind_rows(npath_rows(perf), signed, share) |>
    filter(measure %in% measures) |>
    mutate(run = if_else(.data[[x]] %in% separate, "separate", "grid"),
           measure = factor(measure, measures),
           outcome = factor(outcome_labs[outcome], outcome_labs[outcomes]),
           estimator = factor(estimator, arms))
}

#' measures along n on a log axis in facet strips. several outcomes = measures in rows
#' separate run gets its interval and a dotted join. arms too heavy tailed for the axis left out
fig_npath_strips <- function(s, arms, outcomes = outcome_models,
                             measures = c("RMSE", "|Bias|"),
                             x = s$info$x, block = s$info$specs[1],
                             separate = NULL,
                             breaks = c(250, 1e3, 1e4, 1e5, 1e6),
                             x_lab = expression(italic(n)), title = NULL,
                             legend = "bottom", shade_from = NULL) {
  d <- npath_strip_rows(s, arms, outcomes, measures, x, block, separate) |>
    mutate(xv = .data[[x]])
  grid <- filter(d, run == "grid")
  sep <- filter(d, run == "separate")
  join <- if (nrow(sep)) {
    bind_rows(grid |> group_by(estimator, outcome, measure) |>
                filter(xv == max(xv)) |> ungroup(), sep)
  } else {
    sep
  }
  # free y axes would start where the data start. share runs 0 to 1
  limits <- distinct(d, outcome, measure) |>
    reframe(y = c(0, if (measure == "Share of bias") 1 else 0),
            .by = c(outcome, measure)) |>
    mutate(xv = breaks[1])
  zero <- filter(distinct(d, outcome, measure), measure == "Bias")
  p <- ggplot(d, aes(xv, y, colour = estimator, linetype = estimator,
                     shape = estimator, group = estimator)) +
    geom_blank(data = limits, aes(xv, y), inherit.aes = FALSE)
  if (!is.null(shade_from)) {
    p <- p + annotate("rect", xmin = shade_from, xmax = Inf, ymin = -Inf,
                      ymax = Inf, fill = "#EDE7DB")
  }
  p <- p +
    {if (nrow(zero)) geom_hline(data = zero, aes(yintercept = 0),
                                linewidth = 0.3, colour = "grey40")} +
    geom_ribbon(data = grid, aes(ymin = lo, ymax = hi, fill = estimator),
                colour = NA, alpha = 0.18, show.legend = FALSE) +
    # ggh4x independent scales fail on a layer with no rows
    # so no bar or dotted join w/o a separate run
    {if (nrow(sep)) list(
      geom_linerange(data = sep, aes(ymin = lo, ymax = hi), linewidth = 0.6,
                     alpha = 0.5, show.legend = FALSE),
      geom_line(data = join, linetype = "11", linewidth = 0.4,
                show.legend = FALSE))} +
    geom_line(data = grid, linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 0.8) +
    scale_x_log10(breaks = breaks, labels = format_kM(breaks),
                  expand = expansion(mult = 0.03)) +
    scale_colour_manual(values = npath_strip_pal[arms],
                        labels = arm_labs[arms]) +
    scale_fill_manual(values = npath_strip_pal[arms], guide = "none") +
    scale_linetype_manual(values = npath_strip_lty[arms],
                          labels = arm_labs[arms], guide = "none") +
    scale_shape_manual(values = npath_strip_shape[arms],
                       labels = arm_labs[arms]) +
    labs(x = x_lab, y = NULL, colour = NULL, linetype = NULL, shape = NULL,
         title = title) +
    guides(colour = guide_legend(override.aes = list(size = 1.8))) +
    theme_bw(base_size = 9, base_family = "serif") +
    theme(plot.title = element_text(size = 10, margin = margin(b = 3)),
          plot.title.position = "plot",
          legend.text = element_text(size = 7.5),
          legend.key.width = unit(0.8, "lines"),
          legend.key.height = unit(0.75, "lines"),
          legend.key.spacing.y = unit(0, "pt"),
          legend.margin = margin(2, 4, 2, 4),
          legend.background = element_rect(fill = "white", colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2),
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          strip.placement = "outside",
          # room between panels so last x label of one and first of the next do not run together
          panel.spacing.x = unit(10, "pt"),
          axis.text = element_text(colour = "black"))
  p <- if (is.numeric(legend)) {
    p + theme(legend.position = "inside",
              legend.position.inside = legend,
              legend.justification.inside = c(1, 1))
  } else if (legend == "top") {
    # one row between title and panels flush left w/ the title
    p + theme(legend.position = "top",
              legend.justification.top = "left",
              legend.location = "plot",
              legend.box.spacing = unit(2, "pt"),
              legend.key.spacing.x = unit(6, "pt"))
  } else {
    # one compact row tight under the x axis title
    p + theme(legend.position = "bottom",
              legend.box.spacing = unit(0, "pt"),
              legend.key.spacing.x = unit(6, "pt"))
  }
  p <- if (length(outcomes) > 1) {
    p + ggh4x::facet_grid2(measure ~ outcome, scales = "free_y",
                           independent = "y", switch = "y")
  } else {
    p + facet_wrap(~measure, nrow = 1, scales = "free_y")
  }
  if (!"Share of bias" %in% measures) {
    return(p + ggh4x::facetted_pos_scales(
      y = list(TRUE ~ scale_y_continuous(labels = fmt_ticks))))
  }
  # share read against 0 and .5 and 1 and needs less height than the others
  share_axis <- list(measure == "Share of bias" ~
                       scale_y_continuous(limits = c(0, 1),
                                          breaks = c(0, 0.5, 1),
                                          labels = fmt_ticks),
                     measure != "Share of bias" ~
                       scale_y_continuous(labels = fmt_ticks))
  p <- p + ggh4x::facetted_pos_scales(y = share_axis)
  if (length(outcomes) > 1) {
    p <- p + ggh4x::force_panelsizes(
      rows = ifelse(measures == "Share of bias", 0.7, 1))
  }
  p
}

#' largest MCSE of each measure in fig_npath_strips() for the caption
npath_strips_mcse <- function(s, arms, outcomes = outcome_models,
                              block = s$info$specs[1]) {
  drawn <- filter(s$perf, spec == block, estimator %in% arms,
                  outcome %in% outcomes, !dropped)
  c(rmse = max(drawn$mcse_rmse), bias = max(drawn$mcse_bias),
    share = max(drawn$mcse_share))
}

# sparsity. RMSE as confounders are added (fig-sp-change)

#' geometric mean RMSE per estimator over the sparsity sets by cov structure at each s
#' log gap = log(1 - Rbar/100) so it matches app F tables. checked vs stored balnet vs Hajek
sp_change_rows <- function(s, arms = nc_arms, block = "correct") {
  regime_cells <- list(exact = sparsity_regimes$exact_s$keep,
                       truncated = sparsity_regimes$truncated$keep)
  perf <- filter(s$perf, spec == block, estimator %in% arms)
  if (any(perf$dropped)) stop("a sparsity cell lacks estimates")
  rows <- imap(regime_cells, \(keep, regime) {
    perf[keep(perf), ] |>
      summarise(k = n(), log_rmse = mean(log(rmse)),
                # outcomes share draws so cell errors = perfectly correlated. upper bound on MCSE
                mcse = mean(mcse_rmse / rmse),
                .by = c(covcor, s, estimator)) |>
      mutate(regime = regime)
  }) |>
    list_rbind()

  stored <- s$avg |>
    filter(block == .env$block, outcome == "pooled", theta == "pair",
           A == "balnet0", B == "glmnetcv_hajek", factor %in% "s") |>
    transmute(set, s = as.numeric(value), R_bar)
  gaps <- rows |>
    filter(estimator %in% c("balnet0", "glmnetcv_hajek")) |>
    select(regime, covcor, s, estimator, log_rmse) |>
    pivot_wider(names_from = estimator, values_from = log_rmse) |>
    mutate(set = paste0(if_else(regime == "exact", "exact_s_", "truncated_"),
                        covcor),
           R = 100 * (1 - exp(balnet0 - glmnetcv_hajek))) |>
    inner_join(stored, by = c("set", "s"))
  if (nrow(gaps) != 18 ||
      !isTRUE(all.equal(gaps$R, gaps$R_bar, tolerance = 1e-8))) {
    stop("the sparsity figure no longer matches the stored averages")
  }

  mutate(rows, rmse = exp(log_rmse), lo = exp(log_rmse - 1.96 * mcse),
         hi = exp(log_rmse + 1.96 * mcse))
}

#' arrow chart of RMSE as confounders are added. s = 4 to 64 w/ a circle at 16
fig_sp_change <- function(s, arms = nc_arms, block = "correct") {
  rows <- sp_change_rows(s, arms, block)
  regimes <- c(exact = "Exact balance\np/n = 0.1",
               truncated = "Path truncates\np/n = 0.5, 1")
  structures <- sub("^(.)", "\\U\\1", sparsity_corr, perl = TRUE)
  d <- rows |>
    mutate(regime = factor(regimes[regime], regimes),
           covcor = factor(unname(structures[covcor]), unname(structures)),
           estimator = factor(estimator, rev(arms)))
  ends <- d |>
    select(regime, covcor, estimator, s, rmse) |>
    pivot_wider(names_from = s, values_from = rmse, names_prefix = "s")

  p <- ggplot(d, aes(y = estimator, colour = estimator)) +
    geom_segment(data = ends, aes(x = s4, xend = s64, yend = estimator),
                 arrow = arrow(length = unit(4, "pt"), type = "closed"),
                 linewidth = 0.6) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 1.6, alpha = 0.35) +
    geom_point(data = filter(d, s == 4), aes(x = rmse), shape = 124,
               size = 2.5) +
    geom_point(data = filter(d, s == 16), aes(x = rmse), shape = 21,
               fill = "white", size = 1.4) +
    facet_grid(regime ~ covcor) +
    scale_x_log10(labels = fmt_ticks) +
    scale_y_discrete(labels = arm_labs[arms]) +
    scale_colour_manual(values = npath_strip_pal[arms], guide = "none") +
    labs(x = "RMSE (log scale)", y = NULL,
         title = expression("RMSE as confounders are added, " * italic(s) *
                              " = 4 to 64 (circle at 16)")) +
    theme_bw(base_size = 9, base_family = "serif") +
    theme(plot.title = element_text(size = 10, margin = margin(b = 3)),
          plot.title.position = "plot",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2),
          strip.background = element_blank(),
          strip.text = element_text(size = 8.5),
          axis.text = element_text(colour = "black", size = 7.5))

  k <- distinct(rows, regime, k)
  note <- sprintf(paste(
    "Each point is the geometric mean over %d cells ($n$ = 1000) while",
    "exact balance holds and %d ($n$ = 500 and 1000) once the path",
    "truncates. The largest MCSE of log RMSE is %s."),
    k$k[k$regime == "exact"], k$k[k$regime == "truncated"],
    format(signif(max(rows$mcse), 2)))
  list(plot = p, note = note, height = 3.4)
}

#' after balnet truncates (p/n 0.5 and 1 correct spec). bias and RMSE of the augmented pair
#' and nonzero glmnet PS coefs by cov structure decay and s. means over each row's cells
table_sp_cells <- function(s, s_equi, font = font_serif) {
  perf <- bind_rows(s$full$correct, s_equi$full$correct) |>
    filter(!dropped)
  trunc <- filter(perf, pn >= 0.5)
  bias_arms <- c("balnet0", "abw_balnet0", "glmnetcv_hajek", "aipw_glmnetcv")
  rmse_arms <- c("abw_balnet0", "aipw_glmnetcv")
  keys <- c("covcor", "decay", "s")
  # whatever the corr between a row's cells the mean of their MCSEs bounds the MCSE of the mean
  avg <- trunc |>
    filter(estimator %in% bias_arms) |>
    summarise(bias = mean(bias), rmse = mean(rmse),
              mcse_bias = mean(mcse_bias), mcse_rmse = mean(mcse_rmse),
              k = n(), .by = all_of(c(keys, "estimator")))
  stopifnot(n_distinct(avg$k) == 1)
  # PS fit does not depend on the path so its range spans every p/n of the row
  nnz <- bind_rows(s$diag$correct, s_equi$diag$correct) |>
    filter(estimator == "nnz_ps") |>
    summarise(nnz = paste0(min(centre), "–", max(centre)), .by = all_of(keys))
  wide <- bind_rows(
    transmute(avg, across(all_of(keys)), col = paste0("bias_", estimator),
              val = bias),
    filter(avg, estimator %in% rmse_arms) |>
      transmute(across(all_of(keys)), col = paste0("rmse_", estimator),
                val = rmse)
  ) |>
    pivot_wider(names_from = col, values_from = val) |>
    left_join(nnz, by = keys) |>
    arrange(match(covcor, names(sparsity_corr)), desc(decay), s) |>
    # structure and decay named on their first row only
    mutate(first_cov = !duplicated(covcor),
           first_decay = !duplicated(paste(covcor, decay)),
           covcor = if_else(first_cov, unname(sparsity_corr[covcor]), ""),
           decay = if_else(first_decay, format_plain(decay), ""),
           s = format_plain(s)) |>
    select(-first_cov, -first_decay)

  # cell level facts the text quotes. checked against the stored cells
  exact_bias <- perf |>
    filter(pn == 0.1, estimator == "balnet0") |>
    summarise(m = max(abs(bias)), .by = outcome)
  abw_aipw <- trunc |>
    filter(estimator %in% rmse_arms, covcor %in% c("iid", "ar1")) |>
    select(level, outcome, estimator, rmse) |>
    pivot_wider(names_from = estimator, values_from = rmse)
  stopifnot(all(abw_aipw$abw_balnet0 < abw_aipw$aipw_glmnetcv))
  arm_short <- c(balnet0 = "balnet", abw_balnet0 = "ABW",
                 glmnetcv_hajek = "Hájek", aipw_glmnetcv = "AIPW")
  cols <- c(paste0("bias_", bias_arms), paste0("rmse_", rmse_arms))
  wide <- wide |>
    select(all_of(keys), all_of(cols), nnz) |>
    mutate(across(all_of(cols), \(v) format_signed(v, 3)))
  g <- gt(wide) |>
    cols_label(.list = c(
      list(covcor = "Covariates", decay = md("*ν*"), s = md("*s*"),
           nnz = "Nonzero"),
      set_names(as.list(unname(arm_short[sub("^(bias|rmse)_", "", cols)])),
                cols))) |>
    tab_spanner("Bias", starts_with("bias_")) |>
    tab_spanner("RMSE", starts_with("rmse_")) |>
    tab_spanner("Propensity score", nnz) |>
    cols_align("left", all_of(keys)) |>
    cols_align("right", c(all_of(cols), nnz)) |>
    style_booktabs(font)
  g |>
    tab_source_note(md(paste(
      sprintf(paste("Means over the %d cells of each row ($n$ = 500 and",
                    "1000, $p/n$ = 0.5 and 1, linear and quadratic",
                    "outcomes), at 1000 Monte Carlo replications per cell.",
                    "Largest MCSE: %s (bias), %s (RMSE)."),
              avg$k[1], format_plain(signif(max(avg$mcse_bias), 2)),
              format_plain(signif(max(avg$mcse_rmse), 2))),
      "Nonzero: range of the per-cell median nonzero coefficients of the",
      "glmnet propensity score, over $p/n$ = 0.1 to 1.",
      sprintf(paste("ABW has the lower RMSE in all %d cells with independent",
                    "or AR(1) covariates. While exact balance holds",
                    "($p/n$ = 0.1) the absolute bias of balnet is at most",
                    "%s (linear) and %s (quadratic)."),
              nrow(abw_aipw),
              format_plain(round(exact_bias$m[exact_bias$outcome == "linear"], 3)),
              format_plain(round(exact_bias$m[exact_bias$outcome == "quad1"], 3)))
    )))
}

# p/n grid. a measure vs p/n and avgs by regime (fig-pn-rmse afig-app-pn-bias tbl-pn-adv)

#' share of balnet paths ending above the balance tol per spec and cell (larger of the 2 arms)
pn_trunc_share <- function(app) {
  imap(app$trunc, \(tr, spec) {
    tr |>
      summarise(share = max(k / reps), .by = level) |>
      mutate(spec = spec)
  }) |>
    list_rbind() |>
    inner_join(distinct(app$full[[1]], level, n, pn), by = "level")
}

#' RMSE or bias vs p/n. both specs w/ n = 500 and 1000 side by side (layout of fig_nc())
#' bands +/- 1.96 MCSE (Morris et al. 2019). shaded where half the paths truncate. dotted p = n
pn_lines <- function(s, app, measure = c("rmse", "bias"), arms = nc_arms,
                     shade_label = "balnet truncates", title = NULL) {
  measure <- match.arg(measure)
  mcse_col <- paste0("mcse_", measure)
  perf <- filter(s$perf, estimator %in% arms, !dropped)
  if (!nrow(perf)) stop("no cells of ", paste(arms, collapse = ", "))
  pn_levels <- sort(unique(s$perf$pn))
  n_levels <- sort(unique(s$perf$n))
  # plotmath strips so n is italic as in the text
  spec_lab <- \(sp) factor(if_else(sp == "misspec", "Misspecified", "Correct"),
                           c("Correct", "Misspecified"))
  n_lab <- \(n) factor(paste0("italic(n) == ", n),
                       paste0("italic(n) == ", n_levels))
  d <- perf |>
    mutate(y = .data[[measure]],
           lo = y - 1.96 * .data[[mcse_col]],
           hi = y + 1.96 * .data[[mcse_col]],
           lo = if (measure == "rmse") pmax(lo, 0) else lo,
           spec = spec_lab(spec), n_lab = n_lab(n),
           outcome = factor(outcome_labs[outcome], outcome_labs),
           estimator = factor(estimator, arms))
  # one band per panel from the geometric midpoint below the first level where half truncate
  # to the right edge
  bands <- pn_trunc_share(app) |>
    filter(share >= 0.5) |>
    summarise(from = min(pn), .by = c(spec, n)) |>
    mutate(xmin = sqrt(from * pn_levels[match(from, pn_levels) - 1]),
           xmax = Inf, spec = spec_lab(spec), n_lab = n_lab(n))
  corner <- filter(bands, spec == "Misspecified", n == max(n_levels)) |>
    mutate(outcome = factor(levels(d$outcome)[1], levels(d$outcome)),
           x = xmin, label = shade_label)
  if (is.null(shade_label)) corner <- corner[0, ]
  ggplot(d, aes(pn, y, colour = estimator, shape = estimator,
                linetype = estimator, group = estimator)) +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax,
                                ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "#EDE7DB") +
    geom_vline(xintercept = 1, linetype = "11", linewidth = 0.3,
               colour = "grey55") +
    {if (measure == "bias") geom_hline(yintercept = 0, linewidth = 0.3,
                                       colour = "grey40")} +
    geom_text(data = corner, aes(x = x, y = -Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = -0.7, size = 2.6,
              family = "serif", colour = "grey30") +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = estimator), colour = NA,
                alpha = 0.18, show.legend = FALSE) +
    geom_line(linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 0.9) +
    ggh4x::facet_nested(outcome ~ spec + n_lab, scales = "free_y",
                        labeller = label_parsed) +
    scale_x_log10(breaks = pn_levels, labels = fmt_dot(pn_levels)) +
    scale_y_continuous(labels = fmt_ticks) +
    scale_colour_manual(values = npath_strip_pal[arms],
                        labels = arm_labs[arms]) +
    scale_fill_manual(values = npath_strip_pal[arms], guide = "none") +
    scale_linetype_manual(values = npath_strip_lty[arms], guide = "none") +
    scale_shape_manual(values = npath_strip_shape[arms],
                       labels = arm_labs[arms]) +
    guides(colour = guide_legend(override.aes = list(size = 1.8))) +
    {if (measure == "rmse") expand_limits(y = 0)} +
    labs(x = expression(italic(p) / italic(n)),
         y = if (measure == "rmse") "RMSE" else "Bias",
         colour = NULL, shape = NULL, title = title) +
    theme_bw(base_size = 9, base_family = "serif") +
    theme(plot.title = element_text(size = 10, margin = margin(b = 3)),
          plot.title.position = "plot",
          legend.position = "bottom",
          legend.box.spacing = unit(0, "pt"),
          legend.key.spacing.x = unit(6, "pt"),
          legend.text = element_text(size = 7.5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2),
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          panel.spacing.x = unit(8, "pt"),
          axis.text = element_text(colour = "black"),
          # six levels on a short axis so labels small to keep apart
          axis.text.x = element_text(size = 7))
}

#' largest MCSE of a measure per estimator in pn_lines() for the caption
pn_lines_mcse <- function(s, measure = c("rmse", "bias"), arms = nc_arms) {
  measure <- match.arg(measure)
  drawn <- filter(s$perf, estimator %in% arms, !dropped)
  vapply(arms,
         \(a) max(drawn[[paste0("mcse_", measure)]][drawn$estimator == a]), 0)
}

#' Rbar per pair over all cells and exact and truncated plus sbar of A by regime (tbl-pn-adv)
#' no B share as the text never quotes it and A's share says which side the bias is on
table_pn_regimes <- function(s, pairs, blocks = c("correct", "misspec"),
                             sets = c(all = "All", exact = "Exact balance",
                                      truncated = "Path truncates"),
                             note = NULL, theta = c("pair", "common"),
                             font = font_serif) {
  theta <- match.arg(theta)
  regimes <- setdiff(names(sets), "all")
  avg <- s$avg |>
    filter(set %in% names(sets), block %in% blocks, outcome == "pooled",
           .data$theta == .env$theta) |>
    semi_join(pairs, by = c("A", "B"))
  rows <- avg |>
    select(block, A, B, set, R_bar, s_bar_A) |>
    pivot_wider(names_from = set, values_from = c(R_bar, s_bar_A)) |>
    arrange(match(block, blocks), pair_rank(A, B, pairs)) |>
    # estimator named on its first row only so its pairs read as one
    mutate(est_A = if_else(duplicated(paste(block, A)), "",
                           unname(arm_labs[A])),
           est_B = unname(arm_labs[B]),
           block = unname(block_labs[block])) |>
    select(block, est_A, est_B,
           all_of(paste0("R_bar_", names(sets))),
           all_of(paste0("s_bar_A_", regimes)))
  # cells per spec. a pooled block holds twice as many
  cells <- avg |>
    filter(block != "both") |>
    distinct(block, set, n_cells) |>
    summarise(k = unique(n_cells), .by = set)
  regime_labs <- set_names(as.list(unname(sets[regimes])), regimes)

  gt(rows, groupname_col = "block") |>
    tab_style(cell_text(weight = "bold"), cells_row_groups()) |>
    tab_spanner(md("$\\bar R(\\Theta)$"), starts_with("R_bar_")) |>
    tab_spanner(md("$\\bar s_m(\\Theta)$, $m = A$"),
                starts_with("s_bar_A_")) |>
    cols_label(.list = c(
      list(est_A = "Estimator (A)", est_B = "Comparator (B)"),
      set_names(as.list(unname(sets)), paste0("R_bar_", names(sets))),
      set_names(regime_labs, paste0("s_bar_A_", regimes)))) |>
    fmt_number(starts_with("R_bar_"), decimals = 1) |>
    fmt_number(starts_with("s_bar_"), decimals = 2) |>
    cols_align("left", c(est_A, est_B)) |>
    cols_align("right", -c(est_A, est_B)) |>
    style_booktabs(font) |>
    # hairline above each estimator's first row to split the estimators
    tab_style(cell_borders(sides = "top", color = "#BDBDBD",
                           weight = px(0.5)),
              cells_body(rows = est_A != "" & seq_along(est_A) > 1)) |>
    tab_source_note(md(paste(c(
      sprintf(paste("Cells per specification: %d in all, %d while exact",
                    "balance holds and %d once the path truncates, at",
                    "10 000 Monte Carlo replications."),
              cells$k[cells$set == "all"], cells$k[cells$set == "exact"],
              cells$k[cells$set == "truncated"]),
      note_avg_mcse(avg), note), collapse = " ")))
}

# study 2. interior optima and selector wins (tbl-s2-counts)

# study 2 grids in tbl-grids-s2 order and names
# no folds grid since its ten fold run repeats noise x overlap on fresh draws
s2_grids <- tibble::tribble(
  ~id,                   ~grid,                                   ~factor,
  "cv_overlap",          "Overlap",                               "$c$",
  "cv_noise",            "Outcome noise",                         "$\\kappa$",
  "cv_noise_x_overlap",  "Noise × overlap",                       "$\\kappa$, $c$",
  "cv_alpha_x_overlap",  "Penalty form × overlap",                "$\\alpha$, $c$",
  "cv_n",                "Sample size × overlap × penalty",       "$n$, $c$, $\\alpha$",
  "cv_nbig",             "Large $n$ × overlap",                   "$n$, $c$",
  "cv_p",                "Dimension × overlap × penalty",         "$p$, $c$, $\\alpha$",
  "cv_misspec",          "Misspecification × overlap",            "$c$",
  "cv_sparsity",         "Denser profile",                        "$s$, covariance",
  "cv_sparsity_misspec", "Denser profile under Tan's transformation", "$s$, covariance"
)

s2_selectors <- c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")
s2_est_labs <- c(ipw = "balnet", abw = "ABW")

#' interior optimum = path min over gain below endpoint. win = selector below it by z MCSEs
#' stored interior and beats_floor used other thresholds so both redone. path min must match
s2_cell_flags <- function(s, gain = 0.05, z = 3) {
  key <- c("part", "cell", "outcome", "strength", "est")
  end <- s$perf |>
    filter(arm == "fixed") |>
    select(all_of(key), end = rmse, end_mcse = mcse_rmse, reps)
  path_min <- s$perf |>
    filter(arm == "path") |>
    summarise(path_min = min(rmse), .by = all_of(key))
  selectors <- s$perf |>
    filter(arm %in% s2_selectors) |>
    select(all_of(key), arm, rmse, mcse_rmse) |>
    inner_join(select(end, -reps), by = key) |>
    summarise(win = any(rmse < end - z * sqrt(mcse_rmse^2 + end_mcse^2)),
              lead = arm[which.min(rmse)], .by = all_of(key))
  flags <- end |>
    inner_join(path_min, by = key) |>
    inner_join(selectors, by = key) |>
    mutate(interior = path_min < (1 - gain) * end)
  stored <- inner_join(flags, select(s$cells, all_of(key), stored = path_min),
                       by = key)
  stopifnot(nrow(flags) == nrow(s$cells), nrow(stored) == nrow(flags),
            isTRUE(all.equal(stored$path_min, stored$stored)))
  flags
}

#' kappa as printed in the thesis. 0 1/32 ... 1/2 1
format_kappa <- function(x) {
  stopifnot(all(x %in% c(0, 2^-(0:10))))
  if_else(x %in% c(0, 1), format_plain(x), paste0("1/", 1 / x))
}

#' tbl-s2-counts. share of cells w/ interior optimum and selector win and win given interior
#' whole % since these are counts of cells w/ no MC error. pooled rows weight grids by cells
table_s2_counts <- function(flags, font = font_serif) {
  stopifnot(setequal(flags$id, s2_grids$id))
  shares <- function(d, by) {
    d |>
      # given first since summarise() would otherwise read the overwritten columns
      summarise(given = 100 * mean(win[interior]), cells = n(),
                interior = 100 * mean(interior), win = 100 * mean(win),
                .by = all_of(c(by, "est"))) |>
      pivot_wider(names_from = est, values_from = c(interior, win, given))
  }
  by_grid <- shares(flags, "id") |>
    inner_join(s2_grids, by = "id") |>
    arrange(match(id, s2_grids$id))
  all_grids <- shares(mutate(flags, id = "all"), "id") |>
    mutate(grid = "All grids", factor = "")
  by_kappa <- shares(flags, "strength") |>
    arrange(strength) |>
    mutate(grid = paste("$\\kappa$ =", format_kappa(strength)), factor = "")
  rate_cols <- paste0(rep(c("interior", "win", "given"), each = 2), "_",
                      names(s2_est_labs))
  rows <- bind_rows(A = bind_rows(by_grid, all_grids), B = by_kappa,
                    .id = "panel") |>
    select(panel, grid, factor, cells, all_of(rate_cols))

  outside <- flags |>
    filter(!interior) |>
    summarise(won = sum(win), cells = n(), .by = est)
  outside <- outside[match(names(s2_est_labs), outside$est), ]
  # only the n = 10 000 cells of the large n grid ran fewer reps
  few <- filter(flags, reps < max(reps))
  stopifnot(all(few$id == "cv_nbig"), all(few$part == "n10000"))
  note <- paste(
    "Interior optimum: the path minimum is more than 5% below the",
    "endpoint's RMSE. Selector win: a selector's RMSE is below the",
    "endpoint's by more than three times their combined MCSE; – where a",
    "grid has no interior optimum. Cells per estimator: design cells ×",
    "$\\kappa$ levels × three outcomes; pooled rows weight grids by cells.",
    sprintf(paste("Without an interior optimum a selector won in %d of %d",
                  "cells for balnet and %d of %d for ABW."),
            outside$won[1], outside$cells[1], outside$won[2],
            outside$cells[2]),
    sprintf("%d replications per cell, %d at $n$ = 10 000.",
            max(flags$reps), min(flags$reps))
  )

  gt(rows) |>
    cols_hide(panel) |>
    tab_row_group(md("B. By outcome signal $\\kappa$, all grids"),
                  rows = panel == "B", id = "B") |>
    tab_row_group("A. By grid", rows = panel == "A", id = "A") |>
    row_group_order(c("A", "B")) |>
    tab_style(cell_text(weight = "bold"), cells_row_groups()) |>
    tab_spanner("Interior optimum (%)", starts_with("interior_")) |>
    tab_spanner("Selector win (%)", starts_with("win_")) |>
    tab_spanner("Win given an interior optimum (%)", starts_with("given_")) |>
    cols_label(.list = c(
      list(grid = "Grid", factor = "Factor varied", cells = "Cells"),
      set_names(as.list(rep(unname(s2_est_labs), 3)), rate_cols)
    )) |>
    fmt_markdown(c(grid, factor)) |>
    fmt_number(cells, decimals = 0, sep_mark = "") |>
    fmt_number(all_of(rate_cols), decimals = 0) |>
    sub_missing(all_of(rate_cols), missing_text = "–") |>
    cols_align("left", c(grid, factor)) |>
    cols_align("right", -c(grid, factor)) |>
    style_booktabs(font) |>
    # hairline sets the pooled row apart from the grids it pools
    tab_style(cell_borders(sides = "top", color = "#BDBDBD",
                           weight = px(0.5)),
              cells_body(rows = grid == "All grids")) |>
    tab_source_note(md(paste(note, collapse = " ")))
}

# conclusion. main pairs over every grid (tbl-concl)

#' Rbar of some pairs over every cell of each grid (tbl-concl). each grid's own stored all set
table_rbar_grids <- function(grids, pairs, block = "both", note = NULL,
                             theta = c("pair", "common"), font = font_serif) {
  theta <- match.arg(theta)
  avg <- imap(grids, \(s, grid) {
    s$avg |>
      filter(set == "all", .data$block == .env$block, outcome == "pooled",
             .data$theta == .env$theta) |>
      semi_join(pairs, by = c("A", "B")) |>
      mutate(grid = grid)
  }) |>
    list_rbind()
  wanted <- tidyr::expand_grid(grid = names(grids), pairs)
  absent <- anti_join(wanted, avg, by = c("grid", "A", "B"))
  if (nrow(absent)) {
    stop("no stored average for ",
         paste(unique(paste(absent$grid, absent$A, "vs", absent$B)),
               collapse = "; "))
  }
  pair_cols <- paste0("R_", seq_len(nrow(pairs)))
  pair_labs <- set_names(as.list(pair_name(pairs$A, pairs$B)), pair_cols)
  rows <- avg |>
    mutate(col = pair_cols[match(paste(A, B), paste(pairs$A, pairs$B))]) |>
    select(grid, col, R_bar) |>
    pivot_wider(names_from = col, values_from = R_bar) |>
    # a pair that lost cells avgs over fewer. count shown is the grid's
    # and the note says when a pair skipped some
    left_join(summarise(avg, cells = max(n_cells), .by = grid), by = "grid") |>
    arrange(match(grid, names(grids))) |>
    select(grid, cells, all_of(pair_cols))

  gt(rows) |>
    tab_spanner(md("$\\bar R(\\Theta)$, A vs B"), all_of(pair_cols)) |>
    cols_label(.list = c(list(grid = "Grid", cells = "Cells"), pair_labs)) |>
    fmt_number(all_of(pair_cols), decimals = 1) |>
    fmt_number(cells, decimals = 0, sep_mark = "") |>
    cols_align("left", grid) |>
    cols_align("right", -grid) |>
    style_booktabs(font) |>
    tab_source_note(md(paste(c(
      note_avg_mcse(avg, shares = FALSE), note_dropped(avg, theta), note),
      collapse = " ")))
}
