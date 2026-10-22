# rechecks a built summary from the raw results w/ its own reader and code from the thesis formulas
# not rbar.R or build.R. stored as s$verify and read_summary() warns if an error check failed

library(data.table)

v_outcomes <- c("linear", "quad1", "exp")
v_diag <- "^(lam_|nnz_|smd|ess|wmax|time_|cvloss_|logcvloss_|mse_|msre_|mu_|trunc05|prev|emin|emax|nout)"

#' one row of the verification table. level "error" must pass and "info" is only reported
v_row <- function(check, ok, value = NA_real_, detail = "", level = "error") {
  data.table(check = check, level = level, ok = isTRUE(ok),
             value = as.numeric(value), detail = detail)
}

#' largest abs diff. NA safe and Inf if the NA patterns differ
v_maxdiff <- function(x, y) {
  if (length(x) != length(y)) return(Inf)
  if (!identical(is.na(x), is.na(y))) return(Inf)
  if (all(is.na(x))) return(0)
  max(abs(x - y), na.rm = TRUE)
}

# read (independent of build.R)

#' estimate rows of a section long by outcome. e = tau_hat as tau = 0
v_read <- function(section) {
  parts <- list()
  for (spec_ in names(section$files)) {
    for (file_ in names(section$files[[spec_]])) {
      arms <- section$files[[spec_]][[file_]]
      d <- fread(here::here("results", "aipw", file_),
                 colClasses = list(character = "err"))
      d <- d[estimator %in% arms & !grepl(v_diag, estimator)]
      if (!is.null(section$keep)) {
        rows <- section$keep(d)
        d <- d[rows]
      }
# n = 1M batches are one cell. number reps across them so draws match the stored ones
      if ("batch" %in% names(d)) {
        m_batch <- max(d$sim)
        d[, sim := (batch - 1L) * m_batch + sim]
        d[, batch := NULL]
      }
      oc <- intersect(v_outcomes, names(d))
      ids <- setdiff(names(d), c("err", oc))
      d <- melt(d, id.vars = ids, measure.vars = oc,
                variable.name = "outcome", value.name = "e",
                variable.factor = FALSE)
      d[, `:=`(e = as.numeric(e), spec = spec_)]
      parts[[length(parts) + 1]] <- d
    }
  }
  d <- rbindlist(parts, use.names = TRUE)
  if (!is.null(section$derive)) d <- as.data.table(section$derive(d))
  design <- setdiff(names(d), c("sim", "estimator", "outcome", "e", "spec"))
  list(d = d, design = design)
}

#' balnet path ends per draw and arm. repeated in every outcome col so linear is used
v_read_lam <- function(section) {
  parts <- list()
  for (spec_ in names(section$files)) {
    for (file_ in names(section$files[[spec_]])) {
      d <- fread(here::here("results", "aipw", file_),
                 colClasses = list(character = "err"))
      d <- d[estimator %in% c("lam_end1", "lam_end0")]
      if (!is.null(section$keep)) d <- d[section$keep(d)]
      d[, `:=`(lam = as.numeric(linear), spec = spec_)]
      d[, c("err", intersect(v_outcomes, names(d))) := NULL]
      parts[[length(parts) + 1]] <- d
    }
  }
  d <- rbindlist(parts, use.names = TRUE)
  if (!is.null(section$derive)) d <- as.data.table(section$derive(d))
  d
}

# recompute

#' per cell measures of every arm. all NA if < 2 estimates like arm_measures() in rbar.R
v_perf <- function(d, design) {
  d[, {
    x <- e[is.finite(e)]
    m <- length(x)
    na <- NA_real_
    if (m < 2) {
      list(reps = .N, missing = .N - m, bias = na, empse = na, rmse = na,
           mcse_bias = na, mcse_empse = na, mcse_rmse = na, share = na,
           mcse_share = na)
    } else {
      b <- mean(x)
      M <- mean(x^2)
# share s = b^2 / M. IF = (2b/M)(x - b) - (b^2/M^2)(x^2 - M)
      if_s <- (2 * b / M) * (x - b) - (b^2 / M^2) * (x^2 - M)
      list(reps = .N, missing = .N - m, bias = b, empse = sd(x),
           rmse = sqrt(M), mcse_bias = sd(x) / sqrt(m),
           mcse_empse = sd(x) / sqrt(2 * (m - 1)),
           mcse_rmse = sd(x^2) / (2 * sqrt(M) * sqrt(m)),
           share = b^2 / M, mcse_share = sd(if_s) / sqrt(m))
    }
  }, by = c("spec", design, "estimator", "outcome")]
}

#' paired A and B rows w/ cell missing shares + kept flag + per row influence values
v_pair <- function(d, design, A, B, miss_max) {
  key <- c("spec", design, "sim", "outcome")
  a <- d[estimator == A, c(key, "e"), with = FALSE]
  b <- d[estimator == B, c(key, "e"), with = FALSE]
  run <- intersect(unique(a$spec), unique(b$spec))
  w <- merge(a[spec %in% run], b[spec %in% run], by = key, all = TRUE,
             suffixes = c("_A", "_B"))
  cellkey <- c("spec", design, "outcome")
  w[, `:=`(miss_A = mean(!is.finite(e_A)), miss_B = mean(!is.finite(e_B))),
    by = cellkey]
  w[, kept := miss_A <= miss_max & miss_B <= miss_max]
  cells <- unique(w[, c(cellkey, "miss_A", "miss_B", "kept"), with = FALSE])
  p <- w[kept == TRUE & is.finite(e_A) & is.finite(e_B)]
  p[, `:=`(MA = mean(e_A^2), MB = mean(e_B^2), bA = mean(e_A),
           bB = mean(e_B)), by = cellkey]
# log R_c = (log MA - log MB) / 2. IF(log M) = e^2 / M - 1 and the 1s cancel
  p[, `:=`(i_R = (e_A^2 / MA - e_B^2 / MB) / 2,
           i_sA = (2 * bA / MA) * (e_A - bA) - (bA^2 / MA^2) * (e_A^2 - MA),
           i_sB = (2 * bB / MB) * (e_B - bB) - (bB^2 / MB^2) * (e_B^2 - MB))]
# RMSE_A - RMSE_B paired by draw. delta method IF(RMSE) = e^2 / (2 RMSE)
# its constant drops out of the sd
  pc <- p[, list(reps = .N, rmse_A = sqrt(MA[1]), rmse_B = sqrt(MB[1]),
                 R_c = sqrt(MA[1] / MB[1]),
                 log_R_c_mcse = sd(i_R) / sqrt(.N),
                 rmse_diff = sqrt(MA[1]) - sqrt(MB[1]),
                 rmse_diff_mcse = sd(e_A^2 / (2 * sqrt(MA)) -
                                       e_B^2 / (2 * sqrt(MB))) / sqrt(.N),
                 s_A = bA[1]^2 / MA[1], s_A_mcse = sd(i_sA) / sqrt(.N),
                 s_B = bB[1]^2 / MB[1], s_B_mcse = sd(i_sB) / sqrt(.N)),
          by = cellkey]
  cells <- merge(cells, pc, by = cellkey, all.x = TRUE)
  setkeyv(p, cellkey)                       # fast lookups in v_average()
  list(rows = p, cells = cells)
}

#' Rbar and sbar + MCSEs over the cells `th` from the paired rows
v_average <- function(rows, th, design) {
  k <- nrow(th)
  if (k == 0) return(list(n_cells = 0L, R_bar = NA_real_, R_bar_mcse = NA_real_,
                          s_bar_A = NA_real_, s_bar_A_mcse = NA_real_,
                          s_bar_B = NA_real_, s_bar_B_mcse = NA_real_))
  cellkey <- c("spec", design, "outcome")
  r <- rows[th[, cellkey, with = FALSE], on = cellkey, nomatch = NULL]
# one draw = design level + sim whatever the spec and outcome
  g <- r[, list(gR = sum(i_R), gA = sum(i_sA), gB = sum(i_sB)),
         by = c(design, "sim")]
  v <- g[, list(vR = var(gR) / .N, vA = var(gA) / .N, vB = var(gB) / .N),
         by = design]
  L <- mean(log(th$R_c))
  list(n_cells = k, R_bar = 100 * (1 - exp(L)),
       R_bar_mcse = 100 * exp(L) * sqrt(sum(v$vR)) / k,
       s_bar_A = mean(th$s_A), s_bar_A_mcse = sqrt(sum(v$vA)) / k,
       s_bar_B = mean(th$s_B), s_bar_B_mcse = sqrt(sum(v$vB)) / k)
}

# bootstrap over draws

#' bootstrap sds of Rbar and sbars. draws resampled within level jointly for all cells they feed
v_bootstrap <- function(rows, th, design, B = 200) {
  cellkey <- c("spec", design, "outcome")
  r <- rows[th[, cellkey, with = FALSE], on = cellkey, nomatch = NULL]
  r[, lev := do.call(paste, c(.SD, sep = "|")), .SDcols = design]
  r[, cell := paste(spec, outcome)]
  k <- nrow(th)
  sumL <- numeric(B); sumA <- numeric(B); sumB <- numeric(B)
  for (lv in unique(r$lev)) {
    x <- r[lev == lv]
    sims <- sort(unique(x$sim)); cells <- unique(x$cell)
    at <- cbind(match(x$sim, sims), match(x$cell, cells))
    mk <- function(v) { M <- matrix(0, length(sims), length(cells)); M[at] <- v; M }
    OK <- mk(1); A1 <- mk(x$e_A); A2 <- mk(x$e_A^2); B1 <- mk(x$e_B)
    B2 <- mk(x$e_B^2)
    W <- t(stats::rmultinom(B, length(sims), rep(1, length(sims))))  # B x m
    den <- W %*% OK
    mA <- (W %*% A1) / den; MA <- (W %*% A2) / den
    mB <- (W %*% B1) / den; MB <- (W %*% B2) / den
    sumL <- sumL + rowSums(0.5 * (log(MA) - log(MB)))
    sumA <- sumA + rowSums(mA^2 / MA)
    sumB <- sumB + rowSums(mB^2 / MB)
  }
  c(R_bar = sd(100 * (1 - exp(sumL / k))), s_bar_A = sd(sumA / k),
    s_bar_B = sd(sumB / k))
}

#' groups of cells of a set. whole set or each level of its by factors (values as fmt_value())
v_groups <- function(th0, set) {
  if (is.null(set$by)) {
    return(list(list(f = NA_character_, v = NA_character_,
                     rows = rep(TRUE, nrow(th0)))))
  }
  unlist(lapply(set$by, \(f) lapply(sort(unique(th0[[f]])), \(x)
    list(f = f, v = format(x, scientific = FALSE, drop0trailing = TRUE,
                           trim = TRUE), rows = th0[[f]] == x))),
    recursive = FALSE)
}

#' recompute extras from results alone (path ends + trunc shares + cost of misspec)
#' nc trunc shares (no share col) are not rechecked here
v_extras <- function(s, section, d, design, miss_max, tol) {
  out <- list()
  add <- function(...) out[[length(out) + 1]] <<- v_row(...)
  ex <- s$extra
  if (!is.null(ex$path_ends) || "share" %in% names(ex$trunc)) {
    lam <- v_read_lam(section)
    key <- c("spec", intersect(design, names(lam)))
    w <- dcast(lam, as.formula(paste(paste(c(key, "sim"), collapse = " + "),
                                     "~ estimator")), value.var = "lam")
    if (!is.null(ex$path_ends)) {
      pe <- w[, list(reps = .N,
                     lam1 = mean(lam_end1),
                     lam1_mcse = sd(lam_end1) / sqrt(.N),
                     lam0 = mean(lam_end0),
                     lam0_mcse = sd(lam_end0) / sqrt(.N),
                     lam_sum = mean(lam_end1 + lam_end0),
                     lam_sum_mcse = sd(lam_end1 + lam_end0) / sqrt(.N)),
              by = key]
      st <- as.data.table(ex$path_ends)
      j <- merge(pe, st, by = key, suffixes = c("", ".s"))
      add("path ends: every cell found in both routes",
          nrow(j) == nrow(pe) && nrow(j) == nrow(st), nrow(j),
          sprintf("%d recomputed, %d stored", nrow(pe), nrow(st)))
      vars <- c("reps", "lam1", "lam1_mcse", "lam0", "lam0_mcse", "lam_sum",
                "lam_sum_mcse")
      dmax <- max(sapply(vars, \(v) v_maxdiff(j[[v]], j[[paste0(v, ".s")]])))
      add("path ends: means and MCSEs agree", dmax <= tol, dmax,
          "largest |difference| between the routes")
    }
    if ("share" %in% names(ex$trunc)) {
      st <- as.data.table(ex$trunc)
      k2 <- intersect(key, names(st))
      tr <- w[, list(reps = .N,
                     share = mean(lam_end1 > 1.01e-4 | lam_end0 > 1.01e-4)),
              by = k2]
      tr[, share_mcse := sqrt(share * (1 - share) / reps)]
      j <- merge(tr, st, by = k2, suffixes = c("", ".s"))
      add("truncation: every cell found in both routes",
          nrow(j) == nrow(tr) && nrow(j) == nrow(st), nrow(j),
          sprintf("%d recomputed, %d stored", nrow(tr), nrow(st)))
      dmax <- max(sapply(c("reps", "share", "share_mcse"),
                         \(v) v_maxdiff(j[[v]], j[[paste0(v, ".s")]])))
      add("truncation: shares and MCSEs agree", dmax <= tol, dmax,
          "largest |difference| between the routes")
    }
  }
  if (!is.null(ex$cost) && !is.null(section$cost)) {
    sc <- as.data.table(ex$cost)
    d2 <- d[estimator %in% section$cost$arms]
    d2[, `:=`(estimator = paste0(estimator, "@", spec), spec = "cost")]
    d_cost <- 0; n_rows <- 0; n_found <- 0
    for (cost_arm in section$cost$arms) {
      vp <- v_pair(d2, design, paste0(cost_arm, "@misspec"),
                   paste0(cost_arm, "@correct"), miss_max)
      for (sn in names(section$cost$sets)) {
        set <- section$cost$sets[[sn]]
        th0 <- vp$cells
        if (!is.null(set$keep)) th0 <- th0[set$keep(th0)]
        for (g in v_groups(th0, set)) {
          for (o in c("pooled", intersect(v_outcomes, th0$outcome))) {
            tg <- th0[g$rows]
            if (o != "pooled") tg <- tg[outcome == o]
            a <- v_average(vp$rows, tg[kept == TRUE], design)
            row <- sc[arm == cost_arm & set == sn & outcome == o &
                        (is.na(factor) & is.na(g$f) | factor %in% g$f) &
                        (is.na(value) & is.na(g$v) | value %in% g$v)]
            n_rows <- n_rows + 1
            if (nrow(row) != 1) { d_cost <- Inf; next }
            n_found <- n_found + 1
            if (row$n_cells != a$n_cells ||
                row$n_dropped != nrow(tg) - sum(tg$kept)) d_cost <- Inf
            d_cost <- max(d_cost,
                          v_maxdiff(row$ratio, 1 - a$R_bar / 100),
                          v_maxdiff(row$ratio_mcse, a$R_bar_mcse / 100),
                          v_maxdiff(row$s_bar_A, a$s_bar_A),
                          v_maxdiff(row$s_bar_B, a$s_bar_B))
          }
        }
      }
    }
    add("cost: every stored row recomputed",
        n_found == n_rows && n_rows == nrow(sc), n_found,
        sprintf("%d recomputed, %d stored", n_rows, nrow(sc)))
    add("cost: ratios, MCSEs, shares and cell counts agree", d_cost <= tol,
        d_cost, "largest |difference| between the routes")
  }
  rbindlist(out)
}

# verify one stored summary

#' verify a built summary. B = 0 skips the bootstrap. tol = max diff allowed between routes
verify_summary <- function(s, section, B = 200, tol = 1e-9, seed = 20261007) {
  out <- list()
  add <- function(...) out[[length(out) + 1]] <<- v_row(...)
  miss_max <- s$info$miss_max
  rd <- v_read(section)
  d <- rd$d
  design <- rd$design
  cellkey <- c("spec", design, "outcome")

# per cell perf of every arm
  pv <- v_perf(d, design)
  ps <- as.data.table(s$perf)[, c(cellkey, "estimator", names(pv)[!names(pv) %in%
                                      c(cellkey, "estimator")]), with = FALSE]
  j <- merge(pv, ps, by = c(cellkey, "estimator"), suffixes = c("", ".s"))
  add("perf: every arm-cell found in both routes",
      nrow(j) == nrow(pv) && nrow(j) == nrow(ps), nrow(j),
      sprintf("%d recomputed, %d stored", nrow(pv), nrow(ps)))
  dmax <- max(sapply(c("bias", "empse", "rmse", "mcse_bias", "mcse_empse",
                       "mcse_rmse", "share", "mcse_share"),
                     \(v) v_maxdiff(j[[v]], j[[paste0(v, ".s")]])))
  add("perf: bias, EmpSE, RMSE, share and MCSEs agree", dmax <= tol, dmax,
      "largest |difference| between the routes")
  add("perf: missing counts agree",
      v_maxdiff(j$missing, j$missing.s) == 0, v_maxdiff(j$missing, j$missing.s))
  out <- c(out, list(v_extras(s, section, d, design, miss_max, tol)))

# pairs. per cell R_c and shares then every stored average
  pairs <- as.data.table(s$info$pairs)
  if (nrow(pairs) == 0) {
    add("pairs: none in this section (per-cell measures only)", TRUE, 0,
        level = "info")
    return(rbindlist(out))
  }
  sc <- as.data.table(s$cells)
  sa <- as.data.table(s$avg)
  compared <- unique(c(pairs$A, pairs$B))
  common <- pv[estimator %in% compared,
               list(common = all(missing / reps <= miss_max)), by = cellkey]
  dup <- anyDuplicated(sc[, c("A", "B", "spec", "level", "outcome")])
  add("cells: each cell once per pair", dup == 0, dup)

  blocks <- as.list(s$info$specs)
  names(blocks) <- s$info$specs
  if (length(s$info$specs) > 1) blocks$both <- s$info$specs
  sets <- c(list(all = list()), section$sets)
  d_cells <- 0; d_avg <- 0; n_avg <- 0; n_found <- 0
  for (i in seq_len(nrow(pairs))) {
    pa <- pairs$A[i]; pb <- pairs$B[i]
    vp <- v_pair(d, design, pa, pb, miss_max)
    vc <- merge(vp$cells, common, by = cellkey, all.x = TRUE)
    st <- sc[A == pa & B == pb]
    jc <- merge(vc, st, by = cellkey, suffixes = c("", ".s"))
    if (nrow(jc) != nrow(vc) || nrow(jc) != nrow(st)) d_cells <- Inf
    for (v in c("R_c", "log_R_c_mcse", "s_A", "s_A_mcse", "s_B", "s_B_mcse",
                "rmse_A", "rmse_B", "rmse_diff", "rmse_diff_mcse")) {
      d_cells <- max(d_cells, v_maxdiff(jc[[v]], jc[[paste0(v, ".s")]]))
    }
    d_cells <- max(d_cells, sum(jc$kept != jc$kept.s))
    for (bl in names(blocks)) {
      if (!all(blocks[[bl]] %in% vc$spec)) next
      for (sn in names(sets)) {
        set <- sets[[sn]]
        th0 <- vc[spec %in% blocks[[bl]]]
        if (!is.null(set$keep)) th0 <- th0[set$keep(th0)]
        for (f in set$by) {
          if (anyNA(th0[[f]])) d_avg <- Inf   # a cell in no group of f
        }
        groups <- if (is.null(set$by)) list(list(f = NA_character_,
                                                 v = NA_character_,
                                                 rows = rep(TRUE, nrow(th0))))
        else unlist(lapply(set$by, \(f) lapply(sort(unique(th0[[f]])), \(x)
          list(f = f, v = format(x, scientific = FALSE, drop0trailing = TRUE,
                                 trim = TRUE), rows = th0[[f]] == x))),
          recursive = FALSE)
        for (g in groups) {
          for (th in c("pair", "common")) {
            for (o in c("pooled", intersect(v_outcomes, th0$outcome))) {
              tg <- th0[g$rows]
              if (o != "pooled") tg <- tg[outcome == o]
              use <- tg$kept & (th == "pair" | tg$common %in% TRUE)
              a <- v_average(vp$rows, tg[use], design)
              row <- sa[block == bl & set == sn & theta == th & A == pa &
                          B == pb & outcome == o &
                          (is.na(factor) & is.na(g$f) | factor %in% g$f) &
                          (is.na(value) & is.na(g$v) | value %in% g$v)]
              n_avg <- n_avg + 1
              if (nrow(row) != 1) { d_avg <- Inf; next }
              n_found <- n_found + 1
              if (row$n_cells != a$n_cells ||
                  row$n_dropped != nrow(tg) - sum(use)) d_avg <- Inf
              for (v in c("R_bar", "R_bar_mcse", "s_bar_A", "s_bar_A_mcse",
                          "s_bar_B", "s_bar_B_mcse")) {
                d_avg <- max(d_avg, v_maxdiff(row[[v]], a[[v]]))
              }
            }
          }
        }
      }
    }
# bootstrap of the MCSEs over every kept cell and all outcomes
    if (B > 0) {
      set.seed(seed + i)
      for (bl in names(blocks)) {
        if (!all(blocks[[bl]] %in% vc$spec)) next
        thb <- vc[spec %in% blocks[[bl]] & kept == TRUE]
        if (nrow(thb) == 0) next
        bs <- v_bootstrap(vp$rows, thb, design, B)
        row <- sa[block == bl & set == "all" & theta == "pair" &
                    A == pa & B == pb & outcome == "pooled"]
        ratio <- bs / c(row$R_bar_mcse, row$s_bar_A_mcse, row$s_bar_B_mcse)
        add(sprintf("bootstrap/delta MCSE, %s vs %s, %s", pa, pb, bl),
            all(ratio > 0.8 & ratio < 1.25, na.rm = TRUE), ratio[1],
            sprintf("bootstrap SD / delta MCSE: Rbar %.2f, sbar_A %.2f, sbar_B %.2f (B = %d)",
                    ratio[1], ratio[2], ratio[3], B), level = "info")
      }
    }
  }
  add("cells: R_c, shares, RMSEs, MCSEs and kept flags agree", d_cells <= tol,
      d_cells, "largest |difference| between the routes")
  add("averages: every stored row recomputed",
      n_found == n_avg && n_avg == nrow(sa), n_found,
      sprintf("%d recomputed, %d stored", n_avg, nrow(sa)))
  add("averages: Rbar, sbar, MCSEs, cell counts agree", d_avg <= tol, d_avg,
      "largest |difference| between the routes")

# no cell counted twice
  k <- c("block", "set", "theta", "factor", "value", "A", "B")
  pooled <- sa[outcome == "pooled"]
  per_oc <- sa[outcome != "pooled", list(nc = sum(n_cells), nd = sum(n_dropped)),
               by = k]
  jj <- merge(pooled, per_oc, by = k)
  add("pooled = sum over outcomes (cells and dropped)",
      all(jj$n_cells == jj$nc & jj$n_dropped == jj$nd), nrow(jj))
  if ("both" %in% sa$block) {
    k2 <- c("set", "theta", "factor", "value", "A", "B", "outcome")
    w <- dcast(sa, as.formula(paste(paste(k2, collapse = "+"), "~ block")),
               value.var = "n_cells")
    w <- w[!is.na(both)]
    add("both = correct + misspecified (cells)",
        all(w$both == w$correct + w$misspec), nrow(w))
  }
  part <- sa[!is.na(factor), list(tot = sum(n_cells + n_dropped)),
             by = c("block", "set", "theta", "A", "B", "outcome", "factor")]
  part[, n_fact := .N, by = c("block", "set", "theta", "A", "B", "outcome")]
  part[, same := uniqueN(tot) == 1, by = c("block", "set", "theta", "A", "B",
                                           "outcome")]
  add("by-sets: every factor splits the same cells once", all(part$same),
      nrow(part), "sum over a factor's levels of (cells + dropped) is the same for every factor of a set")
  cm <- sa[theta == "common",
           list(u = uniqueN(n_cells)), by = c("block", "set", "factor",
                                               "value", "outcome")]
  add("common Theta: every pair averaged over the same number of cells",
      all(cm$u == 1), max(cm$u))
  add("Rbar = 100 (1 - exp(mean log R_c))",
      v_maxdiff(sa$R_bar, 100 * (1 - exp(sa$mean_log_R))) <= 1e-9,
      v_maxdiff(sa$R_bar, 100 * (1 - exp(sa$mean_log_R))))
# same arm's sbar over the same cells from different pairs
  sh <- rbind(sa[theta == "common", list(block, set, factor, value, outcome,
                                          arm = A, s = s_bar_A)],
              sa[theta == "common", list(block, set, factor, value, outcome,
                                          arm = B, s = s_bar_B)])
  sh <- sh[, list(spread = max(s) - min(s)),
           by = c("block", "set", "factor", "value", "outcome", "arm")]
  add("sbar of one arm across pairs (common Theta)",
      max(sh$spread, na.rm = TRUE) <= 1e-3, max(sh$spread, na.rm = TRUE),
      "differences come only from pairwise deletion of missing estimates",
      level = "info")
  rbindlist(out)
}

# verify one appendix summary. full tables + diagnostics + approx bias + 800k table from raw files
# MAE MCSE is a bootstrap so only compared to a plain bootstrap as a ratio (info level)

#' one spec of an appendix entry wide (outcomes in cols) w/ its derived cols
v_read_wide <- function(entry, spec, arms = NULL) {
  dir <- if (is.null(entry$dir)) "aipw" else entry$dir
  parts <- lapply(names(entry$files[[spec]]), \(file_) {
    d <- fread(here::here("results", dir, file_),
               colClasses = list(character = "err"))
    keep <- if (is.null(arms)) entry$files[[spec]][[file_]] else arms
    d <- d[estimator %in% keep]
    d[, err := NULL]
    for (o in intersect(v_outcomes, names(d))) {
      set(d, j = o, value = suppressWarnings(as.numeric(d[[o]])))
    }
    d
  })
  d <- rbindlist(parts, use.names = TRUE)
  if (!is.null(entry$derive)) d <- as.data.table(entry$derive(d))
  d
}

#' largest abs diff of cols after joining on `by`. Inf if rows do not match one to one
v_join_diff <- function(x, y, by, cols) {
  x <- as.data.table(x); y <- as.data.table(y)
  j <- merge(x, y, by = by, suffixes = c(".v", ".s"))
  if (nrow(j) != nrow(x) || nrow(j) != nrow(y)) return(Inf)
  max(0, sapply(cols, \(v) v_maxdiff(j[[paste0(v, ".v")]],
                                     j[[paste0(v, ".s")]])))
}

#' verify a built appendix summary. B = plain bootstrap reps for the MAE MCSE check (0 skips)
verify_appendix <- function(a, entry, B = 100, tol = 1e-9, seed = 20261008) {
  out <- list()
  add <- function(...) out[[length(out) + 1]] <<- v_row(...)
  miss_max <- a$info$miss_max

  for (spec_ in names(a$full)) {
    d <- v_read_wide(entry, spec_)
    d <- d[!estimator %in% entry$full$drop]
    oc <- intersect(v_outcomes, names(d))
    design <- setdiff(names(d), c("sim", "estimator", oc))
    long <- melt(d, id.vars = c(design, "sim", "estimator"),
                 measure.vars = oc, variable.name = "outcome",
                 value.name = "e", variable.factor = FALSE)
    pv <- long[, {
      x <- e[is.finite(e)]
      m <- length(x)
      if (m < 2) {
        list(reps = .N, missing = .N - m, bias = NA_real_, empse = NA_real_,
             rmse = NA_real_, mcse_bias = NA_real_, mcse_rmse = NA_real_,
             mcse_empse_k = NA_real_, kurt = NA_real_, mae = NA_real_,
             n_extreme = sum(abs(x) > 10))
      } else {
        dev <- x - mean(x)
        k4 <- mean(dev^4) / mean(dev^2)^2 - 3
        list(reps = .N, missing = .N - m, bias = mean(x), empse = sd(x),
             rmse = sqrt(mean(x^2)), mcse_bias = sd(x) / sqrt(m),
             mcse_rmse = sd(x^2) / (2 * sqrt(mean(x^2)) * sqrt(m)),
             mcse_empse_k = sd(x) / sqrt(2 * (m - 1)) * sqrt(1 + k4 / 2),
             kurt = k4, mae = median(abs(x)), n_extreme = sum(abs(x) > 10))
      }
    }, by = c(design, "estimator", "outcome")]
    st <- as.data.table(a$full[[spec_]])
    dmax <- v_join_diff(pv, st, c(design, "estimator", "outcome"),
                        c("bias", "empse", "rmse", "mcse_bias", "mcse_rmse",
                          "mcse_empse_k", "kurt", "mae", "reps", "missing",
                          "n_extreme"))
    add(sprintf("full (%s): bias, EmpSE, RMSE, MAE, kurtosis, MCSEs, counts agree",
                spec_), dmax <= tol, dmax,
        sprintf("%d arm-cells; largest |difference| between the routes",
                nrow(pv)))
    if (B > 0 && any(!is.na(st$mcse_mae))) {
      set.seed(seed)
      boot <- long[, {
        x <- abs(e[is.finite(e)])
        list(boot_sd = if (length(x) < 2) NA_real_ else
          sd(replicate(B, median(sample(x, replace = TRUE)))))
      }, by = c(design, "estimator", "outcome")]
      j <- merge(boot, st, by = c(design, "estimator", "outcome"))
      ratio <- j$mcse_mae / j$boot_sd
      ratio <- ratio[is.finite(ratio)]
      add(sprintf("full (%s): MAE MCSE against a plain bootstrap", spec_),
          abs(median(ratio) - 1) < 0.1, median(ratio),
          sprintf(paste("stored / plain bootstrap (B = %d): median %.2f,",
                        "middle 90%% %.2f to %.2f"), B, median(ratio),
                  quantile(ratio, 0.05), quantile(ratio, 0.95)),
          level = "info")
    }
  }

  for (spec_ in names(a$diag)) {
    d <- v_read_wide(entry, spec_, arms = c(
      as.vector(outer(as.vector(outer(c("ess", "wmax", "smd"),
                                      c("", "_glm", "_true"), paste0)),
                      c("1", "0"), paste0)),
      "nnz_ps", "nnz_bal1", "nnz_bal0", "nnz_or1", "nnz_or0",
      "lam_end1", "lam_end0", "time_bal", "time_ps", "time_or"))
    oc <- intersect(v_outcomes, names(d))
    design <- setdiff(names(d), c("sim", "estimator", oc))
    long <- melt(d, id.vars = c(design, "sim", "estimator"),
                 measure.vars = oc, variable.name = "outcome",
                 value.name = "v", variable.factor = FALSE)
    per_outcome <- startsWith(long$estimator, "nnz_or")
    long <- long[per_outcome | outcome == oc[1]]
    long[!startsWith(estimator, "nnz_or"), outcome := NA_character_]
    dv <- long[, list(reps = .N, missing = sum(!is.finite(v)),
                      centre = median(v, na.rm = TRUE),
                      lo = quantile(v, 0.25, na.rm = TRUE, names = FALSE),
                      hi = quantile(v, 0.75, na.rm = TRUE, names = FALSE),
                      sd = sd(v, na.rm = TRUE), min = min(v, na.rm = TRUE),
                      max = max(v, na.rm = TRUE)),
               by = c(design, "estimator", "outcome")]
    st <- as.data.table(a$diag[[spec_]])
# NA outcomes do not join in merge() so key them as text
    dv[, outcome := fifelse(is.na(outcome), "-", outcome)]
    st[, outcome := fifelse(is.na(outcome), "-", outcome)]
    dmax <- v_join_diff(dv, st, c(design, "estimator", "outcome"),
                        c("reps", "missing", "centre", "lo", "hi", "sd",
                          "min", "max"))
    add(sprintf("diagnostics (%s): medians, quartiles, SD, range agree", spec_),
        dmax <= tol, dmax,
        sprintf("%d rows; largest |difference| between the routes", nrow(dv)))
    tv <- d[estimator %in% c("lam_end1", "lam_end0"),
            list(k = sum(get(oc[1]) > 1.01e-4, na.rm = TRUE), reps = .N),
            by = c(design, "estimator")]
    ts <- merge(as.data.table(a$trunc[[spec_]]),
                unique(as.data.table(a$diag[[spec_]])[, c("level", design),
                                                      with = FALSE]),
                by = "level")
    dmax <- v_join_diff(tv, ts, c(design, "estimator"), c("k", "reps"))
    add(sprintf("diagnostics (%s): path truncation counts agree", spec_),
        dmax == 0, dmax)
  }

  for (spec_ in names(a$approx)) {
    A <- entry$approx$A; B_ <- entry$approx$B
    d <- v_read_wide(entry, spec_)
    d <- d[estimator %in% c(A, B_)]
    oc <- intersect(v_outcomes, names(d))
    design <- setdiff(names(d), c("sim", "estimator", oc))
    long <- melt(d, id.vars = c(design, "sim", "estimator"),
                 measure.vars = oc, variable.name = "outcome",
                 value.name = "e", variable.factor = FALSE)
    w <- merge(long[estimator == A, c(design, "sim", "outcome", "e"),
                    with = FALSE],
               long[estimator == B_, c(design, "sim", "outcome", "e"),
                    with = FALSE],
               by = c(design, "sim", "outcome"), suffixes = c("_A", "_B"))
    w[, keep := mean(!is.finite(e_A)) <= miss_max &
        mean(!is.finite(e_B)) <= miss_max, by = c(design, "outcome")]
    w <- w[keep == TRUE & is.finite(e_A) & is.finite(e_B)]
    cv <- w[, list(reps = .N, bias_A = mean(e_A), mcse_bias_A = sd(e_A) / sqrt(.N),
                   bias_B = mean(e_B), mcse_bias_B = sd(e_B) / sqrt(.N),
                   bias_diff = mean(e_A - e_B),
                   mcse_bias_diff = sd(e_A - e_B) / sqrt(.N)),
            by = c(design, "outcome")]
    dmax <- v_join_diff(cv, a$approx[[spec_]]$cells, c(design, "outcome"),
                        c("reps", "bias_A", "mcse_bias_A", "bias_B",
                          "mcse_bias_B", "bias_diff", "mcse_bias_diff"))
    add(sprintf("approximation bias (%s): per-cell biases and MCSEs agree",
                spec_), dmax <= tol, dmax)
# averages. centred errors summed over the cells a draw feeds
# then variance over draws within a level and levels added
    avg_of <- function(rows, cells) {
      k <- nrow(cells)
      rows <- copy(rows)
      rows[, `:=`(iA = e_A - mean(e_A), iB = e_B - mean(e_B),
                  iD = (e_A - e_B) - mean(e_A - e_B)),
           by = c(design, "outcome")]
      g <- rows[, list(gA = sum(iA), gB = sum(iB), gD = sum(iD)),
                by = c(design, "sim")]
      v <- g[, list(vA = var(gA) / .N, vB = var(gB) / .N, vD = var(gD) / .N),
             by = design]
      list(n_cells = k, bias_A = mean(cells$bias_A),
           mcse_bias_A = sqrt(sum(v$vA)) / k, bias_B = mean(cells$bias_B),
           mcse_bias_B = sqrt(sum(v$vB)) / k,
           bias_diff = mean(cells$bias_diff),
           mcse_bias_diff = sqrt(sum(v$vD)) / k)
    }
    av <- rbindlist(c(
      list(c(list(outcome = "pooled"), avg_of(w, cv))),
      lapply(intersect(v_outcomes, cv$outcome), \(o) {
        c(list(outcome = o), avg_of(w[outcome == o], cv[outcome == o]))
      })))
    dmax <- v_join_diff(av, a$approx[[spec_]]$avg, "outcome",
                        c("n_cells", "bias_A", "mcse_bias_A", "bias_B",
                          "mcse_bias_B", "bias_diff", "mcse_bias_diff"))
    add(sprintf("approximation bias (%s): averages over cells and MCSEs agree",
                spec_), dmax <= tol, dmax)
  }

  if (!is.null(a$simple)) {
    d <- v_read_wide(entry, "correct")
    oc <- intersect(v_outcomes, names(d))
    design <- setdiff(names(d), c("sim", "estimator", oc))
    long <- melt(d, id.vars = c(design, "sim", "estimator"),
                 measure.vars = oc, variable.name = "outcome",
                 value.name = "e", variable.factor = FALSE)
# every rep counted (none missing)
    sv <- long[, list(m = .N, bias = mean(e), empse = sd(e),
                      rmse = sqrt(mean(e^2)), mcse_bias = sd(e) / sqrt(.N),
                      mcse_empse = sd(e) / sqrt(2 * (.N - 1)),
                      mcse_rmse = sd(e^2) / (2 * sqrt(mean(e^2)) * sqrt(.N))),
               by = c(design, "estimator", "outcome")]
    st <- as.data.table(a$simple)
    for (v in design) if (is.numeric(st[[v]])) set(sv, j = v, value = as.numeric(sv[[v]]))
    dmax <- v_join_diff(sv, st, c(design, "estimator", "outcome"),
                        c("m", "bias", "empse", "rmse", "mcse_bias",
                          "mcse_empse", "mcse_rmse"))
    add("simple: bias, EmpSE, RMSE and MCSEs agree", dmax <= tol, dmax,
        sprintf("%d arm-cells; largest |difference| between the routes",
                nrow(sv)))
  }
  rbindlist(out)
}

# across stored summaries

#' sparsity grid path ends vs the p/n grid where they share the ps (s = 4 nu = 1 iid)
#' correct score once the path truncates (p/n >= 0.5). seed streams differ so draws independent
v_path_end_agreement <- function(sp, pa) {
  sp <- as.data.table(sp)[spec == "correct" & s == 4 & decay == 1 &
                            covcor == "iid" & pn >= 0.5]
  pa <- as.data.table(pa)[spec == "correct"]
  j <- merge(sp, pa, by = c("n", "pn"), suffixes = c("_sp", "_pa"))
  what <- c(treated = "lam1", control = "lam0", sum = "lam_sum")
  rbindlist(lapply(names(what), \(w) {
    v <- what[[w]]
    m_sp <- j[[paste0(v, "_sp")]]
    m_pa <- j[[paste0(v, "_pa")]]
    se <- sqrt(j[[paste0(v, "_mcse_sp")]]^2 + j[[paste0(v, "_mcse_pa")]]^2)
    data.table(n = j$n, pn = j$pn, what = w, sparsity = m_sp, p_axis = m_pa,
               rel_diff = 100 * (m_sp - m_pa) / m_pa, z = (m_sp - m_pa) / se)
  }))
}

#' checks across stored summaries for verification.md. so far sparsity vs p/n path ends
v_across <- function(out_dir) {
  f_sp <- file.path(out_dir, "sparsity.rds")
  f_pa <- file.path(out_dir, "p_axis.rds")
  if (!file.exists(f_sp) || !file.exists(f_pa)) return(data.table())
  sp <- readRDS(f_sp)$extra$path_ends
  pa <- readRDS(f_pa)$extra$path_ends
  if (is.null(sp) || is.null(pa)) return(data.table())
  a <- v_path_end_agreement(sp, pa)
  arm <- a[what != "sum"]
  tot <- a[what == "sum"]
  v_row("path ends, sparsity grid (s = 4, nu = 1, iid) against the p/n grid",
        nrow(tot) > 0 && max(abs(a$z)) < 3, max(abs(a$z)),
        sprintf(paste("correct score, p/n 0.5 and 1 (%d cells): largest",
                      "|difference| %.2f%% (treated or control arm), %.2f%%",
                      "(their sum); largest |z| %.2f (independent draws)"),
                nrow(tot), max(abs(arm$rel_diff)), max(abs(tot$rel_diff)),
                max(abs(a$z))),
        level = "info")
}

# study 2 cv runs (cv.R). second route for build_cv_section() via readr and matrix maths
# sum MCSEs checked by bootstrapping whole reps so it keeps the corr between cells sharing them

#' a cv run file read w/ readr not data.table
v_cv_read <- function(file) {
  path <- here::here("results", "cv", "our dgp", file)
  head <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  chr <- intersect(c("estimator", "err", "covcor"), head)
  types <- do.call(readr::cols,
                   c(setNames(rep(list(readr::col_character()), length(chr)),
                              chr),
                     list(.default = readr::col_double())))
  as.data.table(readr::read_csv(path, col_types = types, progress = FALSE))
}

#' error cols of a cv run file w/ outcome and strength (NA if strength is a design col)
v_cv_cols <- function(d) {
  if (all(v_outcomes %in% names(d))) {
    return(data.table(col = v_outcomes, outcome = v_outcomes,
                      strength = NA_real_))
  }
  col <- grep("^(linear|quad1|exp)_s", names(d), value = TRUE)
  data.table(col = col, outcome = sub("_s.*$", "", col),
             strength = as.numeric(sub("^.*_s", "", col)))
}

#' selector arms of a cv run. the five + cv.balnet ones of each merged file w/ suffix
v_cv_selectors <- function(section) {
  base <- c("cv.bloss", "cv.smd", "cv.inf", "boot.smd", "boot.inf")
  c(base, unlist(lapply(unname(section$merge), \(sfx)
    paste0(c("cv.bloss", "cv.smd", "cv.inf"), sfx))))
}

#' bias EmpSE RMSE share and MCSEs of each col of an error matrix (rows = reps)
v_cv_measures <- function(E) {
  m <- colSums(is.finite(E))
  bias <- colMeans(E, na.rm = TRUE)
  mse <- colMeans(E^2, na.rm = TRUE)
  empse <- apply(E, 2, sd, na.rm = TRUE)
  share <- bias^2 / mse
  inf_share <- sweep(sweep(E, 2, bias) * rep(2 * bias, each = nrow(E)) -
                       sweep(E^2, 2, mse) * rep(share, each = nrow(E)), 2,
                     mse, "/")
  data.table(m = m, bias = bias, empse = empse, rmse = sqrt(mse),
             mcse_bias = empse / sqrt(m),
             mcse_empse = empse / sqrt(2 * (m - 1)),
             mcse_rmse = apply(E^2, 2, sd, na.rm = TRUE) /
               (2 * sqrt(mse) * sqrt(m)),
             share = share,
             mcse_share = apply(inf_share, 2, sd, na.rm = TRUE) / sqrt(m))
}

#' log2 RMSE_A / RMSE_B on same reps. MCSE from IF (e_A^2 / MSE_A - e_B^2 / MSE_B) / 2
v_cv_log2_ratio <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]
  b <- b[ok]
  inf <- (a^2 / mean(a^2) - b^2 / mean(b^2)) / 2
  c(0.5 * log(mean(a^2) / mean(b^2)) / log(2),
    sd(inf) / sqrt(length(a)) / log(2))
}

#' pairs cv.R stores. selectors and grid points vs path end + cv vs boot on smd or inf
v_cv_pairs <- function(selectors) {
  cvc <- grep("^cv\\.(smd|inf)", selectors, value = TRUE)
  rbind(data.table(kind = "selector", A = selectors, B = "fixed"),
        data.table(kind = "path", A = sprintf("path_%02d", 1:29),
                   B = "fixed"),
        data.table(kind = "cv_boot", A = cvc,
                   B = sub("^cv\\.(smd|inf).*$", "boot.\\1", cvc)))
}

#' log RMSE ratio per col w/ paired MCSE + IF values. B NULL gives log RMSE of A
v_cv_log_cols <- function(A, B = NULL, miss_max = 0.01) {
  kept <- colMeans(!is.finite(A)) <= miss_max
  ok <- is.finite(A)
  if (!is.null(B)) {
    kept <- kept & colMeans(!is.finite(B)) <= miss_max
    ok <- ok & is.finite(B)
  }
  ok <- ok & rep(kept, each = nrow(A))
  A[!ok] <- NA
  msa <- colMeans(A^2, na.rm = TRUE)
  inf <- sweep(A^2, 2, msa, "/") - 1
  log_x <- log(msa) / 2
  if (!is.null(B)) {
    B[!ok] <- NA
    msb <- colMeans(B^2, na.rm = TRUE)
    inf <- inf - (sweep(B^2, 2, msb, "/") - 1)
    log_x <- log_x - log(msb) / 2
  }
  inf <- inf / 2
  list(log_x = ifelse(kept, log_x, NA),
       mcse = ifelse(kept, apply(inf, 2, sd, na.rm = TRUE) /
                       sqrt(colSums(ok)), NA),
       inf = inf)
}

#' sum of logs over cols S w/ MCSE from per rep summed IFs (0 if missing) over reps w/ any est
v_cv_set_sum <- function(x, S) {
  S <- S[is.finite(x$log_x[S])]
  present <- rowSums(is.finite(x$inf[, S, drop = FALSE])) > 0
  g <- rowSums(x$inf[present, S, drop = FALSE], na.rm = TRUE)
  c(n = length(S), sum = sum(x$log_x[S]),
    mcse = if (length(S)) sqrt(var(g) / sum(present)) else NA)
}

#' one design cell's stored comparisons recomputed. log ratios per col + sums over sets
v_cv_compare <- function(mats, col, strength, cells, selectors,
                         interior_at = 0.95) {
  pairs <- list()
  sums <- list()
  for (e in c("ipw", "abw")) {
    ce <- cells[est == e]
    ce <- ce[match(col, ce$col)]
    inside <- ce$path_min < interior_at * ce$fixed
    sets <- c(list(list("all", NA_real_, seq_along(col)),
                   list("interior", NA_real_, which(inside)),
                   list("outside", NA_real_, which(!inside))),
              lapply(sort(unique(strength)), \(v)
                list("kappa", v, which(strength == v))))
    add_sums <- function(kind, A, B, x) {
      for (st in sets) {
        r <- v_cv_set_sum(x, st[[3]])
        if (r[["n"]] == 0) next
        sums[[length(sums) + 1]] <<- data.table(
          est = e, kind = kind, A = A, B = B, set = st[[1]], value = st[[2]],
          n_cells = as.integer(r[["n"]]), sum_log = r[["sum"]],
          sum_log_mcse = r[["mcse"]])
      }
    }
    pl <- v_cv_pairs(selectors)
    for (i in seq_len(nrow(pl))) {
      x <- v_cv_log_cols(mats[[paste0(e, "_", pl$A[i])]],
                         mats[[paste0(e, "_", pl$B[i])]])
      pairs[[length(pairs) + 1]] <- data.table(
        col = col, est = e, kind = pl$kind[i], A = pl$A[i], B = pl$B[i],
        log_R = x$log_x, log_R_mcse = x$mcse)[is.finite(log_R)]
      if (pl$kind[i] != "path") add_sums(pl$kind[i], pl$A[i], pl$B[i], x)
    }
    add_sums("arm", "fixed", NA_character_,
             v_cv_log_cols(mats[[paste0(e, "_fixed")]]))
    opt <- mats[[paste0(e, "_fixed")]]
    for (j in seq_along(col)) {
      opt[, j] <- mats[[sprintf("%s_path_%02d", e, ce$k_opt[j])]][, j]
    }
    add_sums("arm", "path_min", NA_character_, v_cv_log_cols(opt))
  }
  list(pairs = rbindlist(pairs), sums = rbindlist(sums))
}

#' every measure of one cell of a cv run recomputed. mats kept for the bootstrap
v_cv_cell <- function(dc, cols, selectors, strength_cell = NA) {
  setorder(dc, estimator, sim)
  rows <- grep("^(ipw|abw)_", unique(dc$estimator), value = TRUE)
  mats <- lapply(setNames(rows, rows), \(r)
    as.matrix(dc[estimator == r, cols$col, with = FALSE]))
  perf <- rbindlist(lapply(rows, \(r)
    cbind(data.table(estimator = r, col = cols$col), v_cv_measures(mats[[r]]))))
  cells <- rbindlist(lapply(c("ipw", "abw"), \(e) {
    rbindlist(lapply(seq_len(nrow(cols)), \(j) {
      rmse_of <- \(r) perf[estimator == r & col == cols$col[j], rmse]
      se_of <- \(r) perf[estimator == r & col == cols$col[j], mcse_rmse]
      end <- paste0(e, "_fixed")
      path <- sprintf("%s_path_%02d", e, 1:30)
      path_rmse <- vapply(path, rmse_of, 0)
      k <- which.min(path_rmse)
      sel <- paste0(e, "_", selectors)
      sel_rmse <- vapply(sel, rmse_of, 0)
      b <- which.min(sel_rmse)
      rb <- v_cv_log2_ratio(mats[[sel[b]]][, j], mats[[end]][, j])
      rm_ <- v_cv_log2_ratio(mats[[path[k]]][, j], mats[[end]][, j])
      fixed <- rmse_of(end)
      data.table(col = cols$col[j], est = e, fixed = fixed, se_f = se_of(end),
                 k_opt = k, path_min = path_rmse[k], se_m = se_of(path[k]),
                 best = selectors[b], best_rmse = sel_rmse[b],
                 se_b = se_of(sel[b]),
                 log2_R_best = rb[1], log2_R_best_mcse = rb[2],
                 log2_R_min = rm_[1], log2_R_min_mcse = rm_[2])
    }))
  }))
  cells[, `:=`(gain_floor = fixed / path_min - 1,
               beats_floor = best_rmse < fixed - 2 * sqrt(se_b^2 + se_f^2),
               sel_vs_fixed = best_rmse / fixed - 1,
               sel_vs_min = best_rmse / path_min - 1)]
  cells[, `:=`(interior = gain_floor > 0.05,
               reaches = beats_floor & sel_vs_min < 0.05)]
  lam <- dc[grepl("^lam[01]_", estimator)]
  picks <- lam[, list(median_lambda = median(get(cols$col[1]), na.rm = TRUE)),
               by = estimator]
  picks[, `:=`(side = ifelse(startsWith(estimator, "lam1"), "treated",
                             "control"),
               arm = sub("^lam[01]_", "", estimator))]
  strength <- ifelse(is.na(cols$strength), strength_cell, cols$strength)
  cmp <- v_cv_compare(mats, cols$col, strength, cells, selectors)
  list(perf = perf, cells = cells, picks = picks, mats = mats,
       pairs = cmp$pairs, sums = cmp$sums)
}

#' check a stored cv summary vs its run files. B = 0 skips bootstrap. n_boot = cells sampled
verify_cv <- function(s, section, B = 200, n_boot = 30, tol = 1e-9,
                      seed = 20261011) {
  out <- list()
  add <- function(...) out[[length(out) + 1]] <<- v_row(...)
  selectors <- v_cv_selectors(section)
  merged <- names(section$merge)
  parts <- setdiff(names(section$files), merged)
  perf <- list(); cells <- list(); picks <- list(); boot <- list()
  pairs <- list(); sums <- list(); boot_sum <- list()
  set.seed(seed)
  for (pt in parts) {
    d <- v_cv_read(section$files[[pt]])
    for (m in merged) {
      dm <- v_cv_read(section$files[[m]])
      dm <- dm[grepl("^(ipw|abw|lam[01])_cv\\.", estimator)]
      dm[, estimator := paste0(estimator, section$merge[[m]])]
      d <- rbind(d, dm)
      rm(dm)
    }
    cols <- v_cv_cols(d)
    for (cl in sort(unique(d$cell))) {
      dc <- d[cell == cl]
      strength_cell <- if ("strength" %in% names(dc)) dc$strength[1] else NA
      r <- v_cv_cell(dc, cols, selectors, strength_cell)
      fix <- function(x) {
        x <- merge(x, cols, by = "col")
        x[is.na(strength), strength := strength_cell]
        x[, `:=`(part = pt, cell = cl)]
      }
      perf[[length(perf) + 1]] <- fix(r$perf)
      cells[[length(cells) + 1]] <- fix(r$cells)
      picks[[length(picks) + 1]] <- r$picks[, `:=`(part = pt, cell = cl)]
      pairs[[length(pairs) + 1]] <- fix(r$pairs)
      sums[[length(sums) + 1]] <- r$sums[, `:=`(part = pt, cell = cl)]
# bootstrap whole reps for the MCSE of a sum over all cols of this cell (one random comparison)
# own seed so the bootstrap below draws the same numbers as before this check was added
      if (B > 0) withr::with_seed(seed + cl, {
        e <- sample(c("ipw", "abw"), 1)
        A <- sample(grep("^(cv|boot)\\.(smd|inf)$", selectors, value = TRUE), 1)
        x <- v_cv_log_cols(r$mats[[paste0(e, "_", A)]],
                           r$mats[[paste0(e, "_fixed")]])
        S <- which(is.finite(x$log_x))
        a <- r$mats[[paste0(e, "_", A)]][, S, drop = FALSE]
        b <- r$mats[[paste0(e, "_fixed")]][, S, drop = FALSE]
        keep <- rowSums(!is.finite(a) | !is.finite(b)) == 0
        a <- a[keep, , drop = FALSE]
        b <- b[keep, , drop = FALSE]
        bs <- replicate(B, {
          i <- sample.int(nrow(a), replace = TRUE)
          sum(log(colMeans(a[i, , drop = FALSE]^2) /
                    colMeans(b[i, , drop = FALSE]^2)) / 2)
        })
        boot_sum[[length(boot_sum) + 1]] <-
          sd(bs) / v_cv_set_sum(x, S)[["mcse"]]
      })
# bootstrap on a few random error cols of this cell
      if (B > 0) {
        for (j in sample(nrow(cols), min(2, nrow(cols)))) {
          for (e in c("ipw", "abw")) {
            rc <- r$cells[est == e & col == cols$col[j]]
            a <- r$mats[[paste0(e, "_", rc$best)]][, j]
            b <- r$mats[[paste0(e, "_fixed")]][, j]
# failed reps left out like the stored measures. by arm for RMSE and by pair for the ratio
            b_end <- b[is.finite(b)]
            ok <- is.finite(a) & is.finite(b)
            a <- a[ok]
            b <- b[ok]
            bs_end <- replicate(B, sqrt(mean(sample(b_end, replace = TRUE)^2)))
            bs_ratio <- replicate(B, {
              i <- sample.int(length(b), replace = TRUE)
              v_cv_log2_ratio(a[i], b[i])[1]
            })
            boot[[length(boot) + 1]] <- data.table(
              rmse_ratio = sd(bs_end) / rc$se_f,
              log2_ratio = sd(bs_ratio) / rc$log2_R_best_mcse)
          }
        }
      }
    }
    rm(d)
  }
  pv <- rbindlist(perf)
  cv <- rbindlist(cells)
  kv <- rbindlist(picks)
  av <- rbindlist(pairs)
  sv <- rbindlist(sums)

  key_p <- c("part", "cell", "estimator", "outcome", "strength")
  ps <- as.data.table(s$perf)
  j <- merge(pv, ps, by = key_p, suffixes = c("", ".s"))
  add("cv perf: every arm-cell found in both routes",
      nrow(j) == nrow(pv) && nrow(j) == nrow(ps), nrow(j),
      sprintf("%d recomputed, %d stored", nrow(pv), nrow(ps)))
  vars <- c("bias", "empse", "rmse", "mcse_bias", "mcse_empse", "mcse_rmse",
            "share", "mcse_share")
  dmax <- max(sapply(vars, \(v) v_maxdiff(j[[v]], j[[paste0(v, ".s")]])))
  add("cv perf: bias, EmpSE, RMSE, share and MCSEs agree", dmax <= tol, dmax,
      sprintf("largest |difference| between the routes %.1e", dmax))
  add("cv perf: replication counts agree", all(j$m == j$m.s),
      sum(j$m != j$m.s), sprintf("%d arm-cells differ", sum(j$m != j$m.s)))

  key_c <- c("part", "cell", "est", "outcome", "strength")
  cs <- as.data.table(s$cells)
  jc <- merge(cv, cs, by = key_c, suffixes = c("", ".s"))
  add("cv cells: every cell found in both routes",
      nrow(jc) == nrow(cv) && nrow(jc) == nrow(cs), nrow(jc),
      sprintf("%d recomputed, %d stored", nrow(cv), nrow(cs)))
  num <- c("fixed", "se_f", "path_min", "se_m", "best_rmse", "se_b",
           "gain_floor", "sel_vs_fixed", "sel_vs_min", "log2_R_best",
           "log2_R_best_mcse", "log2_R_min", "log2_R_min_mcse")
  dmax <- max(sapply(num, \(v) v_maxdiff(jc[[v]], jc[[paste0(v, ".s")]])))
  add("cv cells: RMSEs, MCSEs, ratios and log2 ratios agree", dmax <= tol,
      dmax, sprintf("largest |difference| between the routes %.1e", dmax))
  same <- c("k_opt", "best", "interior", "beats_floor", "reaches")
  n_diff <- sum(sapply(same, \(v) sum(jc[[v]] != jc[[paste0(v, ".s")]])))
  add("cv cells: path minimum, best selector and the three rules agree",
      n_diff == 0, n_diff,
      sprintf("%d cells where the routes choose differently", n_diff))

  ks <- as.data.table(s$picks)
  jk <- merge(kv, ks, by = c("part", "cell", "arm", "side"),
              suffixes = c("", ".s"))
  dk <- v_maxdiff(jk$median_lambda, jk$median_lambda.s)
  add("cv picks: median chosen lambda agrees",
      nrow(jk) == nrow(kv) && nrow(jk) == nrow(ks) && dk <= tol, dk,
      sprintf("%d recomputed, %d stored, largest |difference| %.1e",
              nrow(kv), nrow(ks), dk))

  key_a <- c("part", "cell", "est", "outcome", "strength", "kind", "A", "B")
  as_ <- as.data.table(s$pairs)
  ja <- merge(av, as_, by = key_a, suffixes = c("", ".s"))
  add("cv pairs: every pair in every cell found in both routes",
      nrow(ja) == nrow(av) && nrow(ja) == nrow(as_), nrow(ja),
      sprintf("%d recomputed, %d stored", nrow(av), nrow(as_)))
  da <- max(v_maxdiff(ja$log_R, log(ja$R_c)),
            v_maxdiff(ja$log_R_mcse, ja$log_R_c_mcse))
  add("cv pairs: log R_c and its paired MCSE agree", da <= tol, da,
      sprintf("largest |difference| between the routes %.1e", da))

  key_s <- c("part", "cell", "est", "kind", "A", "B", "set", "value")
  ss <- as.data.table(s$sums)
  js <- merge(sv, ss, by = key_s, suffixes = c("", ".s"))
  add("cv sums: every set found in both routes",
      nrow(js) == nrow(sv) && nrow(js) == nrow(ss) &&
        all(js$n_cells == js$n_cells.s), nrow(js),
      sprintf("%d recomputed, %d stored", nrow(sv), nrow(ss)))
  dsum <- max(v_maxdiff(js$sum_log, js$sum_log.s),
              v_maxdiff(js$sum_log_mcse, js$sum_log_mcse.s))
  add("cv sums: sums of logs and their MCSEs agree", dsum <= tol, dsum,
      sprintf("largest |difference| between the routes %.1e", dsum))
  if (B > 0 && length(boot_sum)) {
    rs <- unlist(boot_sum)
    add("cv bootstrap/delta MCSE (sum over a design cell's columns)",
        min(rs) > 0.75 && max(rs) < 1.33, median(rs),
        sprintf(paste("bootstrap SD / delta MCSE over %d design cells:",
                      "median %.2f, all within %.2f to %.2f (B = %d)"),
                length(rs), median(rs), min(rs), max(rs), B),
        level = "info")
  }

  if (B > 0 && length(boot)) {
    bt <- rbindlist(boot)
    if (nrow(bt) > n_boot) bt <- bt[sample(nrow(bt), n_boot)]
    rng <- range(unlist(bt))
    add("cv bootstrap/delta MCSE (path end RMSE, log2_R_best)",
        rng[1] > 0.75 && rng[2] < 1.33, median(unlist(bt)),
        sprintf(paste("bootstrap SD / delta MCSE over %d sampled cells:",
                      "path end RMSE median %.2f, log2_R_best median %.2f,",
                      "all within %.2f to %.2f (B = %d)"),
                nrow(bt), median(bt$rmse_ratio), median(bt$log2_ratio),
                rng[1], rng[2], B),
        level = "info")
  }
  rbindlist(out)
}

# Wong-Chan entry (cv.R). wc side from the 12 files by the 9 oct draft matrix route not cv.R
# dgp2 side from the twin summary perf table. then matching + draft numbers + MCSE bootstrap

#' paired log RMSE ratios of every col vs col `ref`. MCSE from paired IFs of log MSE
v_wc_paired <- function(err, ref) {
  mse <- colMeans(err^2)
  infl <- sweep(err^2, 2, mse, "/") - 1
  list(log_ratio = 0.5 * (log(mse) - log(mse[ref])),
       mcse = apply((infl - infl[, ref]) / 2, 2, sd) / sqrt(nrow(err)))
}

#' one wc file by the matrix route. per outcome best lambda + best selector vs endpoint + 2 rules
v_wc_file <- function(file, reached, selectors) {
  x <- readRDS(here::here("results", "cv", "wz_replication", file))
  res <- Filter(\(r) is.null(r$err), x$res)
  m <- length(res)
  lam_end <- vapply(res, \(r) as.numeric(r$lam_end), numeric(1))
  share <- colMeans(outer(lam_end, x$lam * (1 + 1e-6), "<="))
  keep <- which(share >= reached)
  k_end <- max(keep)
  satt <- vapply(res, \(r) as.numeric(r$satt), numeric(1))
  mats <- list()
  points <- lapply(c("A", "B"), \(o) {
    truth <- if (o == "A") satt else 0
    e_path <- do.call(rbind, lapply(res, \(r) r$est_path[o, ])) - truth
    e_sel <- do.call(rbind, lapply(res, \(r) r$est_sel[o, selectors])) - truth
    e_end <- e_path[, k_end]
    rmse <- sqrt(colMeans(e_path^2))
    k_opt <- keep[which.min(rmse[keep])]
    at_end <- rmse[k_end] <= rmse[k_opt] * (1 + 1e-12)
    gain <- v_wc_paired(cbind(e_path[, k_opt], e_end), 2)
    sel <- v_wc_paired(cbind(e_sel, e_end), ncol(e_sel) + 1)
    rmse_sel <- sqrt(colMeans(e_sel^2))
    best <- which.min(rmse_sel)
    mcse_rmse <- \(e) sd(e^2) / (2 * sqrt(mean(e^2)) * sqrt(length(e)))
    mcse_end <- mcse_rmse(e_end)
    wins <- rmse_sel < rmse[k_end] -
      3 * sqrt(apply(e_sel, 2, mcse_rmse)^2 + mcse_end^2)
    mats[[o]] <<- list(end = e_end, sel = e_sel)
    data.table(
      outcome = o, sigma = as.numeric(x$cell$sigma), mult = as.numeric(x$cell$c),
      m = m, k_end = k_end, k_opt = if (at_end) k_end else k_opt,
      where = if (at_end) "endpoint" else if (k_opt == 1) "lambda_max" else "interior",
      log_gain = if (at_end) 0 else gain$log_ratio[1],
      mcse_gain = if (at_end) 0 else gain$mcse[1],
      best = selectors[best], log_sel = sel$log_ratio[best],
      mcse_sel = sel$mcse[best],
      interior = rmse[k_opt] < 0.95 * rmse[k_end], win = any(wins),
      rmse_end = rmse[k_end], mcse_end = mcse_end)
  })
  list(points = rbindlist(points), mats = mats)
}

#' dgp2 side from the twin summary perf table. estimates from the arm RMSEs alone
v_wc_dgp2 <- function(d2, overlaps, selectors) {
  p <- as.data.table(d2$perf)[est == "ipw" & outcome == "linear" &
                                overlap %in% overlaps]
  p[, .(
    m = m[arm == "fixed"],
    rmse_end = rmse[arm == "fixed"],
    mcse_end = mcse_rmse[arm == "fixed"],
    k_opt = {
      pk <- .SD[arm == "path"][order(k)]
      i <- which.min(pk$rmse)
      if (pk$rmse[i] >= rmse[arm == "fixed"] * (1 - 1e-12)) 30L else pk$k[i]
    },
    log_gain = {
      pk <- .SD[arm == "path"]
      lg <- log(min(pk$rmse) / rmse[arm == "fixed"])
      if (min(pk$rmse) >= rmse[arm == "fixed"] * (1 - 1e-12)) 0 else lg
    },
    best = {
      ps <- .SD[arm %in% selectors]
      ps <- ps[order(match(arm, selectors))]
      ps$arm[which.min(ps$rmse)]
    },
    log_sel = {
      ps <- .SD[arm %in% selectors]
      log(min(ps$rmse) / rmse[arm == "fixed"])
    },
    interior = min(rmse[arm == "path"]) < 0.95 * rmse[arm == "fixed"],
    win = {
      ps <- .SD[arm %in% selectors]
      any(ps$rmse < rmse[arm == "fixed"] -
            3 * sqrt(ps$mcse_rmse^2 + mcse_rmse[arm == "fixed"]^2))
    }
  ), by = .(overlap, strength)]
}

#' verify the stored Wong-Chan entry
verify_wc <- function(s, section, B = 200, tol = 1e-9, seed = 20261012) {
  out <- list()
  add <- function(...) out[[length(out) + 1]] <<- v_row(...)
  selectors <- s$info$selectors
  stored <- as.data.table(s$points)

# Wong-Chan side
  files <- lapply(section$files, \(f) v_wc_file(f, section$reached, selectors))
  wc <- rbindlist(lapply(files, `[[`, "points"))
  wc[, series := ifelse(outcome == "A", "wc_a", "wc_b")]
  sw <- stored[design != "DGP 2"]
  j <- merge(wc, sw, by.x = c("series", "sigma", "mult"),
             by.y = c("series", "level", "sd_logit"), suffixes = c("", ".s"))
  add("wc points: every Wong-Chan point found in both routes",
      nrow(j) == nrow(wc) && nrow(j) == nrow(sw), nrow(j),
      sprintf("%d recomputed, %d stored", nrow(wc), nrow(sw)))
  num <- c("log_gain", "mcse_gain", "log_sel", "mcse_sel", "rmse_end",
           "mcse_end")
  dmax <- max(sapply(num, \(v) v_maxdiff(j[[v]], j[[paste0(v, ".s")]])))
  add("wc points: log ratios, their paired MCSEs and the endpoint RMSE agree",
      dmax <= tol, dmax,
      sprintf("largest |difference| between the routes %.1e", dmax))
  same <- c("m", "k_end", "k_opt", "where", "best", "interior", "win")
  n_diff <- sum(sapply(same, \(v) sum(j[[v]] != j[[paste0(v, ".s")]])))
  add("wc points: endpoint, best fixed lambda, best selector and the two rules agree",
      n_diff == 0, n_diff,
      sprintf("%d points where the routes choose differently", n_diff))

# dgp2 side
  d2 <- readRDS(here::here("output", "thesis_summaries",
                           paste0(section$dgp2, ".rds")))
  dg <- v_wc_dgp2(d2, section$overlaps, selectors)
  sd2 <- stored[design == "DGP 2"]
  sd2[, overlap := 1 / sd_logit]
  jd <- merge(dg, sd2, by.x = c("overlap", "strength"),
              by.y = c("overlap", "level"), suffixes = c("", ".s"))
  add("wc dgp2: every DGP 2 point found in both routes",
      nrow(jd) == nrow(dg) && nrow(jd) == nrow(sd2), nrow(jd),
      sprintf("%d recomputed, %d stored", nrow(dg), nrow(sd2)))
  dd <- max(v_maxdiff(jd$log_gain, jd$log_gain.s),
            v_maxdiff(jd$log_sel, jd$log_sel.s),
            v_maxdiff(jd$rmse_end, jd$rmse_end.s),
            v_maxdiff(jd$mcse_end, jd$mcse_end.s))
  add("wc dgp2: log ratios from the arms' RMSEs equal the paired ones", dd <= tol,
      dd, sprintf("largest |difference| %.1e", dd))
  n_diff <- sum(jd$k_opt != jd$k_opt.s) + sum(jd$best != jd$best.s) +
    sum(jd$interior != jd$interior.s) + sum(jd$win != jd$win.s) +
    sum(jd$m != jd$m.s)
  add("wc dgp2: best fixed lambda, best selector and the two rules agree",
      n_diff == 0, n_diff,
      sprintf("%d points where the routes choose differently", n_diff))
  twin_cells <- as.data.table(d2$cells)[est == "ipw" & outcome == "linear" &
                                          overlap %in% section$overlaps]
  jc <- merge(sd2, twin_cells, by.x = c("overlap", "level"),
              by.y = c("overlap", "strength"), suffixes = c("", ".c"))
  n_best <- sum(jc$best != jc$best.c)
  add("wc dgp2: best selector as in the twin summary's cells table",
      nrow(jc) == nrow(sd2) && n_best == 0, n_best,
      sprintf("%d points differ", n_best), level = "info")
  add("wc dgp2: the twin summary is the one the entry was built from",
      identical(unname(tools::md5sum(here::here("output", "thesis_summaries",
                                               paste0(section$dgp2, ".rds")))),
                s$info$dgp2$md5),
      detail = paste(section$dgp2, "md5", substr(s$info$dgp2$md5, 1, 8)))

# matching. each wc point vs its dgp2 brackets
  tw <- as.data.table(s$twins)
  agree <- vapply(seq_len(nrow(tw)), \(i) {
    w <- tw[i]
    pool <- sd2[sd2$sd_logit == w$sd_logit]
    tws <- if (w$series == "wc_b") pool[level == 0] else
      pool[level %in% c(w$kappa_lo, w$kappa_hi)]
    (w$wc_interior %in% tws$interior) == w$interior_agree &&
      (w$wc_win %in% tws$win) == w$win_agree &&
      all(tws$r2e >= w$r2e_lo - 1e-12 & tws$r2e <= w$r2e_hi + 1e-12)
  }, logical(1))
  add("wc twins: the bracketing cells and both agreement flags agree",
      all(agree) && nrow(tw) == nrow(sw[series != "wc_b" | level == 1]),
      sum(!agree), sprintf("%d of %d pairs differ", sum(!agree), nrow(tw)))

# numbers the 9 oct draft reported
  draft <- data.table(
    design = c("wc", "wc", "dgp2", "dgp2"),
    series = c("wc_a", "wc_b", "dgp2", "dgp2"),
    level = c(100, 1, 0.125, 0), sd_logit = c(1, 4, 1, 4),
    log_gain = c(-0.008648, -1.405029, -0.005603, -0.911983),
    log_sel = c(-0.010439, -0.993197, -0.006251, -0.895734))
  jk <- merge(draft, stored, by = c("series", "level", "sd_logit"),
              suffixes = c(".d", ""))
  dk <- max(v_maxdiff(jk$log_gain.d, jk$log_gain),
            v_maxdiff(jk$log_sel.d, jk$log_sel))
  add("wc draft: four points as the 9 October draft reported them",
      nrow(jk) == 4 && dk < 1e-5, dk,
      sprintf("largest |difference| to six decimals %.1e", dk), level = "info")

# bootstrap over reps of the best selector log ratio MCSE at 3 wc points
  if (B > 0) {
    set.seed(seed)
    picks <- sample(length(files), 3)
    ratios <- vapply(picks, \(i) {
      o <- sample(c("A", "B"), 1)
      pt <- files[[i]]$points[outcome == o]
      e_sel <- files[[i]]$mats[[o]]$sel[, pt$best]
      e_end <- files[[i]]$mats[[o]]$end
      bs <- replicate(B, {
        r <- sample.int(length(e_end), replace = TRUE)
        0.5 * (log(mean(e_sel[r]^2)) - log(mean(e_end[r]^2)))
      })
      sd(bs) / pt$mcse_sel
    }, numeric(1))
    add("wc bootstrap/delta MCSE (best selector's log ratio)",
        min(ratios) > 0.75 && max(ratios) < 1.33, median(ratios),
        sprintf("bootstrap SD / delta MCSE at 3 points: %s (B = %d)",
                paste(sprintf("%.2f", ratios), collapse = ", "), B),
        level = "info")
  }
  rbindlist(out)
}
