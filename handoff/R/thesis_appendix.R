# appendix full tables from the stored app_<grid> summaries (appendix.R builds them)
# loaded w/ thesis.R by handoff/load.R and uses its helpers. grid labels are in s$info$display

# shared pieces

#' label settings of an appendix summary. defaults n in thousands and other grids plain
app_display <- function(s, spec) {
  d <- s$info$display %||% list()
  by <- d$by
  list(
    x = d$x, by = by,
    x_fmt = switch(d$x_fmt %||% "plain", thousands = format_thousands,
                   plain = format_plain),
    x_title = d$x_title %||% d$x, x_name = d$x_name %||% d$x,
    by_fmt = if (is.null(d$by_fmt)) \(v) paste(by, "=", v) else
      \(v) sprintf(d$by_fmt, v),
    design_labs = d$design_labs,
    spec = unname((d$spec %||% character())[spec] %|NA|% d$spec[1])
  )
}

#' x or y when x is NA or empty
`%|NA|%` <- function(x, y) if (length(x) == 0 || is.na(x)) y else x

#' design cols that split the row groups. ones that vary beyond x and by and the ones kept
#' in stored order so a col fixed by the others (p given n and p/n) does not split again
varying_design <- function(d, design, x, by) {
  kept <- character()
  for (v in setdiff(intersect(design, names(d)), c(x, by))) {
    n_without <- nrow(distinct(d, across(all_of(c(x, by, kept)))))
    if (nrow(distinct(d, across(all_of(c(x, by, kept, v))))) > n_without) {
      kept <- c(kept, v)
    }
  }
  kept
}

#' TRUE unless rendering LaTeX. html markup would print as is there and "&" breaks the table
html_type <- function() !isTRUE(knitr::is_latex_output())

#' thousands separator for a note. nbsp in html and a plain space in LaTeX
thin_space <- function() if (html_type()) "&nbsp;" else " "

#' est w/ MCSE beside it (small type in html). MCSE under half the last digit shows as "<0.001"
estimate_with_mcse <- function(est, se, digits, small = html_type()) {
  lt <- if (small) "&lt;" else "<"
  se_txt <- if_else(se < 0.5 * 10^-digits,
                    paste0(lt, format_plain(10^-digits)),
                    sprintf(paste0("%.", digits, "f"), se))
  wrapped <- if (small) {
    paste0("<span style='font-size:75%'>&nbsp;(", se_txt, ")</span>")
  } else paste0(" (", se_txt, ")")
  paste0(format_signed(est, digits), if_else(is.na(se), "", wrapped))
}

# full performance

#' full perf table for one spec (tbl-app-n-perf-m-* app-pn app-nc)
#' dash = arm misses over 1% of its estimates in that cell
table_full <- function(s, spec = s$info$specs[1], outcome = NULL, keep = NULL,
                       arms = NULL, stack = FALSE, mcse = c("cell", "note"),
                       trunc_note = FALSE, trunc_at = 1.01e-4,
                       measures = c("bias", "empse", "rmse", "mae"),
                       extreme = TRUE, empse_mcse = c("kurtosis", "normal"),
                       x_title = NULL, digits = 3, extreme_at = 10,
                       note = NULL, mcse_small = html_type(),
                       font = font_serif) {
  empse_mcse <- match.arg(empse_mcse)
  mcse <- match.arg(mcse)
  ds <- app_display(s, spec)
  x <- ds$x
  by <- if (stack) NULL else ds$by
  x_title <- x_title %||% ds$x_title
  miss_max <- s$info$miss_max
  perf <- s$full[[spec]]
  if (!is.null(keep)) perf <- perf[keep(perf), ]
  arms <- intersect(names(arm_labs),
                    intersect(arms %||% perf$estimator, perf$estimator))
  p <- filter(perf, estimator %in% arms,
              outcome %in% (.env$outcome %||% outcome_models))
  measure_labs <- c(bias = "Bias", empse = "EmpSE", rmse = "RMSE",
                    mae = "MAE")
  se_col <- c(bias = "mcse_bias", rmse = "mcse_rmse", mae = "mcse_mae",
              empse = if (empse_mcse == "kurtosis") "mcse_empse_k" else
                "mcse_empse")
  measures <- intersect(names(measure_labs), measures)

  long <- map(measures, \(k) {
    shown <- if (mcse == "cell") {
      estimate_with_mcse(p[[k]], p[[se_col[[k]]]], digits, mcse_small)
    } else format_signed(p[[k]], digits)
    mutate(p, Measure = measure_labs[[k]],
           val = if_else(dropped, "—", shown))
  }) |>
    list_rbind()
  extreme_arms <- unique(p$estimator[p$n_extreme > 0])
  if (extreme && length(extreme_arms)) {
    long <- bind_rows(long, p |>
                        filter(estimator %in% extreme_arms) |>
                        mutate(Measure = "Extreme",
                               val = as.character(n_extreme)))
  }

  rows <- varying_design(p, s$info$design, x, by)
  design_label <- \(v) coalesce(unname((ds$design_labs %||% character())[v]),
                                v)
  group_label <- \(d) {
    lbl <- unname(outcome_labs[d$outcome])
    for (v in rows) lbl <- paste0(lbl, " · ", design_label(v), " = ",
                                  format_plain(d[[v]]))
    lbl
  }
  levels_x <- p |>
    distinct(across(all_of(c(by, x)))) |>
    arrange(across(all_of(c(by, x))))
  column_id <- \(d) if (is.null(by)) paste0("x_", d[[x]]) else
    paste0("b_", d[[by]], "|x_", d[[x]])
  cols <- column_id(levels_x)
  measure_order <- c(unname(measure_labs[measures]), "Extreme")
  wide <- long |>
    mutate(col = column_id(long),
           oc = factor(outcome, intersect(outcome_models, outcome)),
           Group = group_label(long),
           Arm = factor(unname(arm_labs[estimator]), unname(arm_labs[arms])),
           Measure = factor(Measure, measure_order)) |>
    arrange(oc, across(all_of(rows)), Arm, Measure) |>
    mutate(Group = factor(Group, unique(Group))) |>
    select(Group, Arm, Measure, col, val) |>
    pivot_wider(names_from = col, values_from = val) |>
    arrange(Group, Arm, Measure) |>
    # arm named on its first measure row only
    mutate(Arm = if_else(Measure == measure_order[1], as.character(Arm), ""),
           Measure = as.character(Measure)) |>
    select(Group, Arm, Measure, all_of(cols))

  g <- if (n_distinct(wide$Group) > 1) {
    gt(wide, groupname_col = "Group")
  } else {
    gt(select(wide, -Group))
  }
  g <- g |>
    cols_label(.list = c(list(Arm = "", Measure = ""),
                         set_names(as.list(ds$x_fmt(levels_x[[x]])), cols))) |>
    cols_align("left", c(Arm, Measure)) |>
    cols_align("right", all_of(cols))
  if (is.null(by)) {
    g <- tab_spanner(g, md(x_title), all_of(cols))
  } else {
    for (v in unique(levels_x[[by]])) {
      g <- tab_spanner(g, ds$by_fmt(v), starts_with(paste0("b_", v, "|")))
    }
    g <- cols_label(g, Measure = md(x_title))
  }
  # cells are html already. identity stops gt escaping them
  if (mcse_small) g <- fmt(g, all_of(cols), fns = identity)

  reps <- unique(p$reps)
  of_reps <- if (length(reps) == 1) {
    paste0(", of ", format(reps, big.mark = thin_space(), scientific = FALSE))
  } else ""
  missing_txt <- if (all(p$missing == 0)) "No estimate is missing." else
    note_missing(p, miss_max)
  largest_txt <- if (mcse == "note") {
    kept <- filter(p, !dropped)
    paste0("Largest MCSE: ", paste(map_chr(measures, \(k) sprintf(
      "%s (%s)", format_plain(signif(max(kept[[se_col[[k]]]]), 2)),
      measure_labs[[k]])), collapse = ", "), ".")
  }
  trunc_txt <- if (trunc_note) {
    cell_keys <- c(by, x, rows)
    note_truncation(s$trunc[[spec]], distinct(p, level, across(all_of(
      cell_keys))), cell_keys, c(by, ds$x_name, map_chr(rows, design_label)),
      trunc_at)
  }
  notes <- c(
    if ("mae" %in% measures) "MAE: median of |τ̂ − τ|.",
    largest_txt,
    missing_txt,
    if (extreme && length(extreme_arms)) sprintf(
      "Extreme: estimates with |τ̂ − τ| > %s%s, kept in every measure.",
      format_plain(extreme_at), of_reps),
    paste0("τ = ", format_plain(s$info$tau), "."),
    trunc_txt,
    note
  )
  g |>
    style_booktabs(font) |>
    tab_options(data_row.padding = px(1)) |>
    tab_style(cell_text(whitespace = "nowrap"), cells_body()) |>
    sub_missing(everything(), missing_text = "") |>
    tab_source_note(md(paste(notes, collapse = " ")))
}

#' compact table for many-factor grids (sparsity grid app-sp). one row per cell
#' estimates only. largest MCSE per measure goes in the note
table_full_cells <- function(s, spec = s$info$specs[1], outcome,
                             columns, groups = c("n", "pn"),
                             stub = c("covcor", "decay", "s"),
                             arm_short = c(balnet0 = "balnet",
                                           abw_balnet0 = "ABW",
                                           glmnetcv_hajek = "Hájek",
                                           aipw_glmnetcv = "AIPW"),
                             ps_nnz = FALSE, digits = 3, note = NULL,
                             font = font_serif) {
  ds <- app_display(s, spec)
  measure_labs <- c(bias = "Bias", rmse = "RMSE")
  se_col <- c(bias = "mcse_bias", rmse = "mcse_rmse")
  p <- filter(s$full[[spec]], .data$outcome %in% .env$outcome)
  by_outcome <- length(outcome) > 1
  col_ids <- unlist(imap(columns, \(arms, k) paste0(k, "_", arms)),
                    use.names = FALSE)
  col_arm <- set_names(unlist(columns, use.names = FALSE), col_ids)
  if (by_outcome) {
    col_ids <- unlist(map(outcome, \(o) paste0(o, "_", names(col_arm))))
    col_arm <- set_names(rep(unname(col_arm), length(outcome)), col_ids)
  }
  long <- imap(columns, \(arms, k) {
    p |>
      filter(estimator %in% arms) |>
      mutate(col = paste0(k, "_", estimator),
             val = if_else(dropped, "—", format_signed(.data[[k]], digits)))
  }) |>
    list_rbind()
  if (by_outcome) long <- mutate(long, col = paste0(.data$outcome, "_", col))
  wide <- long |>
    select(level, all_of(c(groups, stub)), col, val) |>
    pivot_wider(names_from = col, values_from = val)
  if (ps_nnz) {
    nnz <- filter(s$diag[[spec]], estimator == "nnz_ps")
    stopifnot(!anyDuplicated(nnz$level))
    wide <- left_join(wide, transmute(nnz, level,
                                      nnz_ps = format_plain(centre)),
                      by = "level")
  }
  design_label <- \(v) coalesce(unname((ds$design_labs %||% character())[v]),
                                v)
  corr_order <- names(sparsity_corr)
  wide <- arrange(wide, across(all_of(groups)), match(covcor, corr_order),
                  desc(decay), s)
  group_lab <- do.call(paste, c(map(groups, \(v) {
    paste(design_label(v), "=", format_plain(wide[[v]]))
  }), sep = ", "))
  wide <- wide |>
    mutate(Group = factor(group_lab, unique(group_lab)),
    # covcor and decay shown on their first row only
    decay_shown = if_else(duplicated(paste(Group, covcor, decay)), "",
                          format_plain(decay)),
    covcor = if_else(duplicated(paste(Group, covcor)), "",
                     unname(sparsity_corr[covcor])),
    decay = decay_shown, s = format_plain(s)) |>
    select(Group, all_of(stub), all_of(col_ids), any_of("nnz_ps"))

  g <- gt(wide, groupname_col = "Group") |>
    cols_label(.list = c(
      list(covcor = "Covariates", decay = md("*ν*"), s = md("*s*")),
      set_names(as.list(unname(arm_short[col_arm[col_ids]])), col_ids),
      if (ps_nnz) list(nnz_ps = "Nonzero")
    )) |>
    cols_align("left", all_of(stub)) |>
    cols_align("right", c(all_of(col_ids), any_of("nnz_ps")))
  for (o in if (by_outcome) outcome else "") {
    for (k in names(columns)) {
      ids <- paste0(if (by_outcome) paste0(o, "_"), k, "_", columns[[k]])
      lab <- if (by_outcome) {
        paste0(measure_labs[[k]], ", ", tolower(outcome_labs[[o]]))
      } else {
        measure_labs[[k]]
      }
      g <- tab_spanner(g, lab, all_of(ids))
    }
  }
  if (ps_nnz) g <- tab_spanner(g, "Propensity score", nnz_ps)

  shown <- filter(p, estimator %in% unlist(columns), !dropped)
  largest <- paste(map_chr(names(columns), \(k) sprintf(
    "%s (%s)", format_plain(signif(max(shown[[se_col[[k]]]]), 2)),
    measure_labs[[k]])), collapse = ", ")
  shorts <- arm_short[unique(unlist(columns))]
  long_names <- shorts != arm_labs[names(shorts)]
  notes <- c(
    if (any(long_names)) paste0(paste(sprintf(
      "%s: %s", shorts[long_names], arm_labs[names(shorts)[long_names]]),
      collapse = "; "), "."),
    paste0("Largest MCSE: ", largest, "."),
    if (all(p$missing == 0)) "No estimate is missing." else
      note_missing(p, s$info$miss_max),
    if (ps_nnz) paste(
      "Propensity score: median nonzero coefficients at lambda.min, one fit",
      "for both arms and all outcomes."),
    paste0("τ = ", format_plain(s$info$tau), "."),
    note
  )
  g |>
    style_booktabs(font) |>
    tab_options(data_row.padding = px(1)) |>
    tab_style(cell_text(whitespace = "nowrap"), cells_body()) |>
    sub_missing(everything(), missing_text = "") |>
    tab_source_note(md(paste(notes, collapse = " ")))
}

#' note on balnet paths ending above the tolerance in the shown cells
#' over 4 cells truncated gives a count. else each cell's treated and control counts
note_truncation <- function(trunc, cells, keys, key_names,
                            trunc_at = 1.01e-4) {
  if (is.null(trunc) || !nrow(trunc)) return(NULL)
  trunc <- filter(trunc, level %in% cells$level)
  sci <- \(v) sub("e([-+])0*", "e\\1", format(v, scientific = TRUE))
  big <- \(v) format(v, big.mark = thin_space(), scientific = FALSE,
                     trim = TRUE)
  located <- \(dd) {
    do.call(paste, c(map2(keys, key_names, \(col, n) {
      paste(n, "=", format_plain(dd[[col]]))
    }), sep = ", "))
  }
  tr <- trunc |>
    inner_join(cells, by = "level") |>
    pivot_wider(names_from = estimator, values_from = k) |>
    filter(lam_end1 > 0 | lam_end0 > 0) |>
    arrange(across(all_of(keys)))
  if (!nrow(tr)) {
    sprintf("No balnet path ends above λ = %s.", sci(trunc_at))
  } else if (nrow(tr) > 4) {
    sprintf(paste("balnet paths ending above λ = %s: in %d of %d cells, up",
                  "to %s (treated) and %s (control) per cell."),
            sci(trunc_at), nrow(tr), n_distinct(trunc$level),
            big(max(tr$lam_end1)), big(max(tr$lam_end0)))
  } else {
    sprintf(paste("balnet paths ending above λ = %s, treated and control:",
                  "%s; none in the other cells."),
            sci(trunc_at),
            paste(sprintf("%s: %s and %s", located(tr), big(tr$lam_end1),
                          big(tr$lam_end0)), collapse = "; "))
  }
}

# fit diagnostics

#' fit diagnostics for one spec. median (quartiles) per cell (tbl-app-n-diag-m-* app-pn app-sp)
#' balnet nnz goes in the note when its median is the same in every cell
table_diagnostics <- function(s, spec = s$info$specs[1], outcome = NULL,
                              keep = NULL, stack = FALSE, times = TRUE,
                              sets = c("bal", "glm", "true"),
                              fits = c("ps", "or"),
                              families = c("nnz", "ess", "wmax", "smd"),
                              digits = c(nnz = 0, ess = 0, wmax = 1, smd = 3),
                              nnz_bal = c("note", "row"), trunc_at = 1.01e-4,
                              note = NULL, small = html_type(),
                              font = font_serif) {
  nnz_bal <- match.arg(nnz_bal)
  ds <- app_display(s, spec)
  x <- ds$x
  by <- if (stack) NULL else ds$by
  set_labs <- c(bal = "balnet", glm = "glmnet", true = "Oracle")
  d <- s$diag[[spec]]
  trunc <- s$trunc[[spec]]
  if (!is.null(keep)) {
    d <- d[keep(d), ]
    # trunc note only counts the shown cells
    trunc <- filter(trunc, level %in% d$level)
  }
  machine <- s$session[[spec]]
  # balance logged on the basis the fits use. under misspec that is not the true X
  basis <- if (grepl("^Misspecified", ds$spec %||% "")) " (X̃)"
  oc <- outcome %||% intersect(outcome_models, d$outcome)
  d <- filter(d, is.na(outcome) | outcome %in% oc)
  times_d <- filter(d, startsWith(estimator, "time_"))
  arm_labs_tc <- c("1" = "Treated", "0" = "Control")
  family_labs <- c(nnz = "Nonzero coefficients", ess = "ESS",
                   wmax = "Largest weight", smd = "Max |SMD|")

  # balnet nnz as one note line when its median is constant
  bal_nnz <- filter(d, estimator %in% c("nnz_bal1", "nnz_bal0"))
  bal_row <- nnz_bal == "row" ||
    (nrow(bal_nnz) > 0 && n_distinct(bal_nnz$centre) > 1)

  # row key = family fit arm (+ outcome for the outcome models)
  d <- d |>
    mutate(
      fam = case_when(startsWith(estimator, "nnz") ~ "nnz",
                      startsWith(estimator, "ess") ~ "ess",
                      startsWith(estimator, "wmax") ~ "wmax",
                      startsWith(estimator, "smd") ~ "smd"),
      set = case_when(estimator == "nnz_ps" ~ "ps",
                      startsWith(estimator, "nnz_or") ~ "or",
                      startsWith(estimator, "nnz_bal") ~ "bal",
                      grepl("_glm", estimator) ~ "glm",
                      grepl("_true", estimator) ~ "true",
                      .default = "bal"),
      arm = if_else(estimator == "nnz_ps", "both",
                    sub(".*([01])$", "\\1", estimator))
    ) |>
    filter(!is.na(fam), fam %in% families,
           set %in% c(sets, fits),
           !(set == "bal" & fam == "nnz" & !bal_row))
  fit_label <- \(fam, set, outcome) case_when(
    set == "ps" ~ "glmnet PS",
    set == "or" & length(oc) > 1 ~ paste0("Lasso OR (",
                                          tolower(outcome_labs[outcome]), ")"),
    set == "or" ~ "Lasso OR",
    fam == "nnz" & set == "bal" ~ "balnet",
    .default = unname(set_labs[set])
  )
  fit_order <- c("glmnet PS", "balnet",
                 paste0("Lasso OR (", tolower(outcome_labs[oc]), ")"),
                 "Lasso OR", unname(set_labs[sets]))

  # cell = centre at the family digits w/ quartiles in small type
  format_digits <- \(v, k) unname(mapply(\(vv, kk) {
    formatC(vv, format = "f", digits = kk, big.mark = " ")
  }, v, rep_len(k, length(v))))
  cell <- \(centre, lo, hi, sd, stat, k) {
    spread <- if_else(stat == "median",
                      paste0(format_digits(lo, k), "–", format_digits(hi, k)),
                      format_digits(sd, k))
    shown <- if (small) {
      paste0(format_digits(centre, k), "<span style='font-size:75%'>&nbsp;(",
             spread, ")</span>")
    } else paste0(format_digits(centre, k), " (", spread, ")")
    # dash if missing in every rep (glmnet weights when its scores hit 0 or 1) like perf tables
    if_else(is.na(centre), "—", shown)
  }
  rows <- varying_design(d, s$info$design, x, by)
  design_label <- \(v) coalesce(unname((ds$design_labs %||% character())[v]),
                                v)
  levels_x <- d |>
    distinct(across(all_of(c(by, x)))) |>
    arrange(across(all_of(c(by, x))))
  column_id <- \(dd) if (is.null(by)) paste0("x_", dd[[x]]) else
    paste0("b_", dd[[by]], "|x_", dd[[x]])
  cols <- column_id(levels_x)
  group_label <- \(dd) {
    lbl <- unname(family_labs[dd$fam])
    for (v in rows) lbl <- paste0(lbl, " · ", design_label(v), " = ",
                                  format_plain(dd[[v]]))
    lbl
  }
  wide <- d |>
    mutate(col = column_id(d),
           Fit = fit_label(fam, set, outcome),
           Arm = unname(c(arm_labs_tc, both = "Both")[arm]),
           val = cell(centre, lo, hi, sd, stat, unname(digits[fam])),
           Group = group_label(d),
           fam = factor(fam, families),
           Fit = factor(Fit, intersect(fit_order, Fit)),
           Arm = factor(Arm, c("Both", "Treated", "Control"))) |>
    arrange(fam, across(all_of(rows)), Fit, Arm) |>
    mutate(Group = factor(Group, unique(Group))) |>
    select(Group, Fit, Arm, col, val) |>
    pivot_wider(names_from = col, values_from = val) |>
    arrange(Group, Fit, Arm) |>
    mutate(Fit = if_else(!duplicated(paste(Group, Fit)),
                         as.character(Fit), ""),
           Arm = as.character(Arm)) |>
    select(Group, Fit, Arm, all_of(cols))

  g <- gt(wide, groupname_col = "Group") |>
    cols_label(.list = c(list(Fit = "", Arm = ""),
                         set_names(as.list(ds$x_fmt(levels_x[[x]])), cols))) |>
    cols_align("left", c(Fit, Arm)) |>
    cols_align("right", all_of(cols))
  if (is.null(by)) {
    g <- tab_spanner(g, md(ds$x_title), all_of(cols))
  } else {
    for (v in unique(levels_x[[by]])) {
      g <- tab_spanner(g, ds$by_fmt(v), starts_with(paste0("b_", v, "|")))
    }
    g <- cols_label(g, Arm = md(ds$x_title))
  }
  if (small) g <- fmt(g, all_of(cols), fns = identity)

  stat <- d$stat[1]
  reps <- unique(d$reps)
  of_reps <- if (length(reps) == 1) {
    paste0(" over the ", format(reps, big.mark = thin_space(),
                                scientific = FALSE), " replications")
  } else " over the replications"
  bal_txt <- if (!bal_row && nrow(bal_nnz)) sprintf(paste(
    "balnet keeps a median of %s covariates per arm in every cell",
    "(range %s–%s)."), format_plain(bal_nnz$centre[1]),
    format_plain(min(bal_nnz$min)), format_plain(max(bal_nnz$max)))
  # path ends are about balnet so only tables w/ balnet rows mention them
  trunc_txt <- if ("bal" %in% sets) {
    note_truncation(trunc, distinct(d, level, across(all_of(c(by, x, rows)))),
                    c(by, x, rows),
                    c(by, ds$x_name, map_chr(rows, design_label)), trunc_at)
  }
  time_txt <- if (times && nrow(times_d)) {
    time_range <- \(e) {
      v <- times_d$centre[times_d$estimator == e]
      if (!length(v)) return(NA)
      paste(unique(trimws(formatC(signif(range(v), 2), format = "fg",
                                  digits = 2))), collapse = "–")
    }
    parts <- c(time_bal = "balnet (both arms)",
               time_ps = "glmnet propensity score (with CV)",
               time_or = paste("lasso outcome models (both arms, every",
                               "outcome, with CV)"))
    parts <- parts[names(parts) %in% times_d$estimator]
    paste0(
      if (stat == "median") "Median" else "Mean",
      " wall-clock time per fit, range over the cells (s): ",
      paste(sprintf("%s %s", parts, map_chr(names(parts), time_range)),
            collapse = "; "),
      ". Timed inside parallel workers with single-threaded BLAS",
      if (!is.null(machine)) paste0(" on ", paste(machine, collapse = ", ")),
      "; timings depend on the machine and load."
    )
  }
  n_missing <- sum(d$missing)
  # captions give median and quartiles so the note only says what differs
  notes <- c(
    if (stat != "median") paste0("Mean (SD)", of_reps, "."),
    "No MCSE, as these describe the fits.",
    if ("nnz" %in% families && length(fits)) paste0(
      "Nonzero coefficients at lambda.min: ",
      paste(c(ps = "the glmnet propensity score, one fit for both arms",
              or = "the lasso outcome model per arm")[fits],
            collapse = " and "), "."),
    if ("nnz" %in% families) bal_txt,
    if ("ess" %in% families) "ESS: Kish's, per arm.",
    if ("wmax" %in% families) paste(
      "Largest weight: the arm's weights scaled to sum to n."),
    if ("smd" %in% families) paste0(
      "Max |SMD|: over the covariates the fits use", basis,
      ", weighted arm mean against the full-sample mean in full-sample SDs."),
    if (any(c("ess", "wmax", "smd") %in% families) && "true" %in% sets) {
      "Oracle: weights from the true propensity score."
    },
    if (n_missing > 0) sprintf("%d diagnostic values missing.", n_missing),
    trunc_txt,
    time_txt,
    note
  )
  g |>
    style_booktabs(font) |>
    tab_options(data_row.padding = px(1)) |>
    tab_style(cell_text(whitespace = "nowrap"), cells_body()) |>
    sub_missing(everything(), missing_text = "") |>
    tab_source_note(md(paste(notes, collapse = " ")))
}

# path truncation on the n x c grid

#' per cent of reps whose balnet path truncates (an arm ends above 1.01e-4) per n x c cell
#' MCSE at most 1.6 points at m = 1000 so whole per cents (tbl-app-nc-trunc)
table_nc_trunc <- function(s, digits = 0, note = NULL, font = font_serif) {
  m <- s$info$reps
  tr <- s$extra$trunc |>
    mutate(Spec = factor(block_labs[spec],
                         block_labs[intersect(names(block_labs), spec)]),
           col = paste0("x_", format_plain(overlap)),
           val = formatC(100 * trunc, format = "f", digits = digits))
  levels_c <- sort(unique(tr$overlap))
  cols <- paste0("x_", format_plain(levels_c))
  rows <- tr |>
    select(Spec, n, col, val) |>
    pivot_wider(names_from = col, values_from = val) |>
    arrange(Spec, n) |>
    mutate(n = format_plain(n)) |>
    select(Spec, n, all_of(cols))
  mcse <- max(sqrt(tr$trunc * (1 - tr$trunc) / m))
  gt(rows, groupname_col = "Spec", rowname_col = "n") |>
    cols_label(.list = set_names(as.list(format_plain(levels_c)), cols)) |>
    tab_spanner(md("Overlap *c*"), all_of(cols)) |>
    tab_stubhead(md("*n*")) |>
    cols_align("right", all_of(cols)) |>
    style_booktabs(font) |>
    tab_style(cell_text(align = "left", whitespace = "nowrap"),
              cells_stub()) |>
    tab_source_note(md(paste(c(
      sprintf(paste("Truncated: the path ends above λ = 1.01e-4. Largest",
                    "MCSE: %.1f percentage points."), 100 * mcse),
      note
    ), collapse = " ")))
}

# Wong-Chan vs DGP 2 for study 2. balnet only. stored entry wc_compare from cv.R sec 6
# matched on the true ps logit SD and R2_e (share of control outcome var the true ps explains)

wc_series_labs <- c(dgp2 = "DGP 2, misspecified grid",
                    wc_a = "Wong–Chan, model A",
                    wc_b = "Wong–Chan, model B")
# DGP 2 in the balancing blue of the line figs and W-C in the glmnet red
# model B is a different outcome so its point is open and not joined
wc_pal <- c(dgp2 = "#2166AC", wc_a = "#B22222", wc_b = "#B22222")

#' points drawn. model B only at sigma = 1 (paper design. larger sigma has R2_e under 1e-4)
wc_compared <- function(s) {
  filter(s$points, series != "wc_b" | level == 1)
}

#' largest MCSE of a log ratio per design over the drawn points (for captions)
wc_compare_mcse <- function(s) {
  d <- wc_compared(s) |>
    mutate(mcse = pmax(mcse_gain, mcse_sel))
  c(wc = max(d$mcse[d$design != "DGP 2"]),
    dgp2 = max(d$mcse[d$design == "DGP 2"]))
}

#' W-C vs DGP 2 fig (afig-app-wc). RMSE ratios to the endpoint vs R2_e. one column per logit SD
#' pseudo-log x so kappa = 0 is on the axis. bars 95% MC intervals (Morris et al. 2019 sec 6.2)
fig_wc_compare <- function(s, title = NULL) {
  pts <- wc_compared(s)
  column_lab <- c(
    `1` = "Logit SD 1 (Wong–Chan ×1, DGP 2 c = 1)",
    `2` = "Logit SD 2 (Wong–Chan ×2, DGP 2 c = 0.5)",
    `4` = "Logit SD 4 (Wong–Chan ×4, DGP 2 c = 0.25)")
  measure_lab <- c(gain = "Best fixed λ ÷ endpoint",
                   sel = "Best selector ÷ endpoint")
  long <- bind_rows(
    transmute(pts, series, sd_logit, r2e, measure = "gain", lr = log_gain,
              mcse = mcse_gain),
    transmute(pts, series, sd_logit, r2e, measure = "sel", lr = log_sel,
              mcse = mcse_sel)) |>
    mutate(ratio = exp(lr), lo = exp(lr - 1.96 * mcse),
           hi = exp(lr + 1.96 * mcse),
           series = factor(series, names(wc_series_labs)),
           column = factor(column_lab[as.character(sd_logit)], column_lab),
           measure = factor(measure_lab[measure], measure_lab))
  x_tr <- scales::transform_pseudo_log(sigma = 1e-4, base = 10)
  ggplot(long, aes(r2e, ratio, colour = series)) +
    geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey40") +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.4,
                   show.legend = FALSE) +
    geom_line(data = filter(long, series != "wc_b"), linewidth = 0.5,
              show.legend = FALSE) +
    geom_point(aes(shape = series, fill = series), size = 1.7, stroke = 0.6) +
    facet_grid(measure ~ column, scales = "free_y", switch = "y") +
    scale_x_continuous(transform = x_tr, breaks = c(0, 1e-3, 1e-2, 0.1, 1),
                       labels = c("0", "0.001", "0.01", "0.1", "1")) +
    scale_y_continuous(transform = "log2", breaks = c(1 / 4, 1 / 2, 1, 2),
                       labels = c("1/4", "1/2", "1", "2")) +
    scale_colour_manual(values = wc_pal, labels = wc_series_labs) +
    scale_fill_manual(values = c(dgp2 = "#2166AC", wc_a = "#B22222",
                                 wc_b = "white"),
                      labels = wc_series_labs) +
    scale_shape_manual(values = c(dgp2 = 21, wc_a = 22, wc_b = 22),
                       labels = wc_series_labs) +
    guides(colour = guide_legend(override.aes = list(size = 1.8))) +
    labs(x = expression(paste("Confounding strength ", italic(R)[e]^2,
                              " (pseudo-log scale)")),
         y = NULL, colour = NULL, fill = NULL, shape = NULL, title = title) +
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
          strip.text = element_text(size = 8.5),
          strip.placement = "outside",
          panel.spacing.x = unit(8, "pt"),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(size = 7))
}

#' kappa as a fraction
wc_kappa_lab <- function(x) {
  vapply(x, \(v) if (v %in% c(0, 1)) format_plain(v) else
    paste0("1/", round(1 / v)), "")
}

#' R2_e to 2 sig digits. model B (0.00027) to 1
wc_r2e_lab <- function(x) {
  vapply(x, \(v) {
    if (v == 0) return("0")
    formatC(v, format = "fg", digits = if (v < 0.001) 1 else 2, flag = "#")
  }, "")
}

#' ratio to 2 decimals or the range over two bracketing DGP 2 cells
wc_ratio_lab <- function(lo, hi) {
  a <- sprintf("%.2f", exp(lo))
  b <- sprintf("%.2f", exp(hi))
  if_else(a == b, a, paste0(a, "–", b))
}

#' rule verdict yes or no or the two verdicts of bracketing cells
wc_flag_lab <- function(lo, hi) {
  yn <- \(v) if_else(v, "yes", "no")
  if_else(lo == hi, yn(lo), paste0(yn(lo), "–", yn(hi)))
}

#' matched pairs table (atbl-app-wc). one row per W-C point grouped by logit SD
#' do best fixed lambda and a selector beat the endpoint in both and do the study 2 rules agree
table_wc_twins <- function(s, font = font_serif) {
  tw <- s$twins
  group_lab <- c(
    `1` = "Propensity logit SD 1 (Wong–Chan ×1, DGP 2 c = 1)",
    `2` = "Propensity logit SD 2 (Wong–Chan ×2, DGP 2 c = 0.5)",
    `4` = "Propensity logit SD 4 (Wong–Chan ×4, DGP 2 c = 0.25)")
  span <- \(lo, hi, lab) if_else(lo == hi, lab(lo),
                                 paste0(lab(lo), "–", lab(hi)))
  rows <- tw |>
    mutate(group = factor(group_lab[as.character(sd_logit)], group_lab),
           wc = paste0(if_else(series == "wc_b", "B", "A"), ", $\\sigma$ = ",
                       sigma, " (", wc_r2e_lab(r2e), ")"),
           d2 = paste0("$\\kappa$ = ", span(kappa_lo, kappa_hi, wc_kappa_lab),
                       " (", span(r2e_lo, r2e_hi, wc_r2e_lab), ")"),
           gain_wc = wc_ratio_lab(wc_gain, wc_gain),
           gain_d2 = wc_ratio_lab(d2_gain_lo, d2_gain_hi),
           sel_wc = wc_ratio_lab(wc_sel, wc_sel),
           sel_d2 = wc_ratio_lab(d2_sel_lo, d2_sel_hi),
           int_wc = wc_flag_lab(wc_interior, wc_interior),
           int_d2 = wc_flag_lab(d2_interior_lo, d2_interior_hi),
           win_wc = wc_flag_lab(wc_win, wc_win),
           win_d2 = wc_flag_lab(d2_win_lo, d2_win_hi)) |>
    select(group, wc, d2, gain_wc, gain_d2, sel_wc, sel_d2, int_wc, int_d2,
           win_wc, win_d2)
  off <- filter(tw, !interior_agree | !win_agree)
  exceptions <- if (nrow(off)) {
    off |>
      mutate(txt = paste0(
        "logit SD ", sd_logit, ", $\\sigma$ = ", sigma, " (",
        if_else(!interior_agree,
                paste0("interior optimum ", wc_flag_lab(wc_interior, wc_interior),
                       " against ", wc_flag_lab(d2_interior_lo, d2_interior_hi)),
                ""),
        if_else(!interior_agree & !win_agree, "; ", ""),
        if_else(!win_agree,
                paste0("selector win ", wc_flag_lab(wc_win, wc_win),
                       " against ", wc_flag_lab(d2_win_lo, d2_win_hi)),
                ""),
        ")")) |>
      pull(txt) |>
      paste(collapse = " and ")
  }
  mc <- wc_compare_mcse(s)
  note <- paste(c(
    paste("Ratios: RMSE at the best fixed $\\lambda$ on the path, and of the",
          "best of the five shared selectors, over RMSE at the path",
          "endpoint; below 1 the endpoint is beaten. Each Wong–Chan cell",
          "is set against the cells of DGP 2's misspecified grid at the",
          "same logit SD whose $R^2_e$ (in brackets) bracket its own, model",
          "B against $\\kappa$ = 0; a DGP 2 entry spans its two cells.",
          "Interior optimum and selector win as defined in Performance",
          "measures, applied to both designs."),
    sprintf(paste("The designs agree on the interior optimum in %d of %d",
                  "matched cells and on the selector win in %d of %d%s."),
            sum(tw$interior_agree), nrow(tw), sum(tw$win_agree), nrow(tw),
            if (nrow(off)) paste0("; the exception: ", exceptions) else ""),
    sprintf(paste("Largest MCSE of a log ratio %.3f (Wong–Chan),",
                  "%.3f (DGP 2)."), mc[["wc"]], mc[["dgp2"]])
  ), collapse = " ")
  gt(rows, groupname_col = "group") |>
    tab_style(cell_text(weight = "bold"), cells_row_groups()) |>
    tab_spanner(md("Matched cells ($R^2_e$)"), c(wc, d2)) |>
    tab_spanner(md("Best fixed $\\lambda$"), c(gain_wc, gain_d2)) |>
    tab_spanner("Best selector", c(sel_wc, sel_d2)) |>
    tab_spanner("Interior optimum", c(int_wc, int_d2)) |>
    tab_spanner("Selector win", c(win_wc, win_d2)) |>
    cols_label(wc = "Wong–Chan", d2 = "DGP 2", gain_wc = "W–C",
               gain_d2 = "DGP 2", sel_wc = "W–C", sel_d2 = "DGP 2",
               int_wc = "W–C", int_d2 = "DGP 2", win_wc = "W–C",
               win_d2 = "DGP 2") |>
    fmt_markdown(c(wc, d2)) |>
    cols_align("left", c(wc, d2)) |>
    cols_align("right", c(gain_wc, gain_d2, sel_wc, sel_d2)) |>
    cols_align("center", c(int_wc, int_d2, win_wc, win_d2)) |>
    style_booktabs(font) |>
    tab_source_note(md(note))
}
