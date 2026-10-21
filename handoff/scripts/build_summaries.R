# builds checks and stores the summaries thesis.R reads. open the .Rproj then source this file
# then build_summaries() for all or eg build_summaries("overlap"). B = 0 skips the bootstraps

source(here::here("handoff", "R", "rbar.R"))
source(here::here("handoff", "R", "sections.R"))
source(here::here("handoff", "R", "build.R"))
source(here::here("handoff", "R", "cv.R"))
source(here::here("handoff", "R", "appendix.R"))
source(here::here("handoff", "R", "verify.R"))

#' build verify and store some sections then rewrite rbar_all.csv and verification.md
build_summaries <- function(ids = c(names(sections), names(appendix_sections),
                                    names(cv_sections), names(wc_sections)),
                            B = 200) {
  known <- c(names(sections), names(appendix_sections), names(cv_sections),
             names(wc_sections))
  unknown <- setdiff(ids, known)
  if (length(unknown)) {
    stop("not a section in sections.R or cv.R: ",
         paste(unknown, collapse = ", "))
  }
  out_dir <- here::here("output", "thesis_summaries")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# build verify and store each section
# 1 to 2 GB per section. ~6 min for the 8 sections + ~5 for the appendix entries (windows R 4.5)
  for (id in ids) {
    message("Building ", id, " ...")
    started <- Sys.time()
    appendix <- id %in% names(appendix_sections)
    cv <- id %in% names(cv_sections)
    wc <- id %in% names(wc_sections)
    s <- if (appendix) {
      build_appendix(appendix_sections[[id]], id = id)
    } else if (cv) {
      build_cv_section(cv_sections[[id]], id = id)
    } else if (wc) {
      build_wc_section(wc_sections[[id]], id = id)
    } else {
      build_section(sections[[id]], id = id)
    }
    if (!all(s$checks$ok)) {
      print(s$checks)
      stop(id, ": a check failed, nothing stored for this section")
    }
    message("  verifying ...")
    s$verify <- if (appendix) {
      verify_appendix(s, appendix_sections[[id]], B = min(B, 100))
    } else if (cv) {
      verify_cv(s, cv_sections[[id]], B = B)
    } else if (wc) {
      verify_wc(s, wc_sections[[id]], B = B)
    } else {
      verify_summary(s, sections[[id]], B = B)
    }
    s$info$verified <- all(s$verify$ok[s$verify$level == "error"])
    saveRDS(s, file.path(out_dir, paste0(id, ".rds")))
# nc (no pairs) and appendix and cv entries have no averages
    if (NROW(s$avg)) {
      readr::write_csv(s$avg, file.path(out_dir, paste0(id, "_averages.csv")))
    }
    message(sprintf("  done in %s; %s",
                    format(round(Sys.time() - started, 1)),
                    if (s$info$verified) "verified" else "VERIFICATION FAILED"))
    rm(s)
    gc(verbose = FALSE)
  }

# all averages in one file + the report
# read back from disk so it covers every section stored not only this run
  stored <- intersect(known,
                      sub("\\.rds$", "", list.files(out_dir, "\\.rds$")))
  all_avg <- list()
  report <- c("# Verification of the stored summaries", "",
              paste0("Report written ", format(Sys.time(), "%Y-%m-%d %H:%M"),
                     " (", R.version.string, "). Each section lists the ",
                     "pairing checks run before storing, then the checks of ",
                     "verify.R. [x] passed; (info) never stops a build."), "")
  for (id in stored) {
    s <- readRDS(file.path(out_dir, paste0(id, ".rds")))
    if (NROW(s$avg)) all_avg[[id]] <- cbind(section = id, s$avg)
    r <- rbind(
      data.frame(check = s$checks$check, level = "error", ok = s$checks$ok,
                 detail = s$checks$detail),
      as.data.frame(s$verify)[c("check", "level", "ok", "detail")])
    bad <- !r$ok & r$level == "error"
    report <- c(report,
                paste0("## ", id, if (any(bad)) " (FAILED)" else " (passed)"),
                "", paste0("Built ", format(s$info$built, "%Y-%m-%d %H:%M"),
                           "; files: ",
                           paste(s$info$files$file, collapse = ", "), "."),
                "",
                sprintf("- [%s] %s%s: %s", ifelse(r$ok, "x", " "),
                        ifelse(r$level == "info", "(info) ", ""), r$check,
                        r$detail), "")
  }
# checks across two summaries (sparsity vs p/n path ends where they share a ps)
  across <- v_across(out_dir)
  if (nrow(across)) {
    report <- c(report, "## Across sections", "",
                sprintf("- [%s] %s%s: %s", ifelse(across$ok, "x", " "),
                        ifelse(across$level == "info", "(info) ", ""),
                        across$check, across$detail), "")
  }
  readr::write_csv(do.call(rbind, all_avg),
                   file.path(out_dir, "rbar_all.csv"))
  report_path <- file.path(out_dir, "verification.md")
  writeLines(report, report_path)
  message("Report: ", report_path)
  invisible(report_path)
}
