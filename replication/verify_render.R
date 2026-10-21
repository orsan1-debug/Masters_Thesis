# renders the thesis from the stored summaries in a temp copy so nothing here is overwritten
#   Rscript replication/verify_render.R   a few min. needs Quarto and XeLaTeX. reads no results/

root  <- here::here()
copy  <- file.path(tempdir(), "thesis_copy")
qmd   <- file.path("thesis", "AAAAA Working Thesis.qmd")
files <- c(qmd, file.path("thesis", "references.bib"), file.path("handoff", "load.R"),
           file.path("handoff", "R", list.files(file.path(root, "handoff", "R"), "\\.R$")),
           file.path("R", "dgp.R"),
           file.path("output", "thesis_summaries",
                     list.files(file.path(root, "output", "thesis_summaries"))))
for (d in unique(dirname(files))) dir.create(file.path(copy, d), recursive = TRUE, showWarnings = FALSE)
stopifnot(all(file.copy(file.path(root, files), file.path(copy, files), copy.date = TRUE)))
invisible(file.create(file.path(copy, ".here")))   # the copy is the root for here::here()

quarto <- Sys.which("quarto")
if (!nzchar(quarto)) {
  quarto <- Filter(file.exists, c("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto",
                                  "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
                                  "C:/Program Files/Quarto/bin/quarto.exe"))[1]
}
if (is.na(quarto)) stop("Quarto not found")

started <- Sys.time()
log <- withr::with_dir(copy, withr::with_envvar(c(QUARTO_R = R.home("bin")),
  suppressWarnings(system2(quarto, c("render", shQuote(qmd), "--to", "pdf", "-M", "keep-tex:true"),
                           stdout = TRUE, stderr = TRUE))))
out <- here::here("replication", "verification")
dir.create(file.path(out, "render_tex"), recursive = TRUE, showWarnings = FALSE)
writeLines(log, file.path(out, "render_log.txt"))
# keep the tex and figures in render_tex/ so the PDF can be built elsewhere w/ xelatex
tex <- list.files(file.path(copy, "thesis"), "\\.tex$|_files$", full.names = TRUE)
invisible(file.copy(tex, file.path(out, "render_tex"), recursive = TRUE, overwrite = TRUE))
pdf <- file.path(copy, "thesis", c("AAAAA-Working-Thesis.pdf", "AAAAA Working Thesis.pdf"))
pdf <- pdf[file.exists(pdf)][1]
if (is.na(pdf)) stop("no PDF, see replication/verification/render_log.txt")
invisible(file.copy(pdf, file.path(out, "thesis_render.pdf"), overwrite = TRUE, copy.mode = FALSE))

# verification/compare_render.py then compares this PDF w/ the stored one
reference <- here::here("output", "rendered", "thesis", "AAAAA-Working-Thesis.pdf")
cat("Quarto", system2(quarto, "--version", stdout = TRUE),
    "\nrendered in", format(Sys.time() - started, digits = 2),
    "\nmd5 render   ", unname(tools::md5sum(file.path(out, "thesis_render.pdf"))),
    "\nmd5 reference", unname(tools::md5sum(reference)), "\n")
