# compares the verify_render.R pdf w/ the stored one. display numbers + all numbers + page words
#   Rscript replication/verification/compare_render.R [stored.pdf] [render.pdf]   needs pdftotext

out    <- Sys.getenv("OUT_DIR", here::here("replication", "verification"))
args   <- commandArgs(TRUE)
stored <- if (length(args) > 0) args[1] else here::here("output", "rendered", "thesis", "AAAAA-Working-Thesis.pdf")
render <- if (length(args) > 1) args[2] else file.path(out, "thesis_render.pdf")
if (!nzchar(Sys.which("pdftotext"))) stop("pdftotext not found. it comes w/ poppler")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

split_words <- function(x) if (nzchar(x)) strsplit(x, " ", fixed = TRUE)[[1]] else character()

#' lines of each page. NFKC and dashes to hyphens. page numbers dropped later as breaks can differ
pdf_pages <- function(pdf, layout) {
  txt <- system2("pdftotext", c(if (layout) "-layout", "-enc", "UTF-8", shQuote(pdf), "-"), stdout = TRUE)
  if (!is.null(attr(txt, "status"))) stop("pdftotext failed on ", pdf)
  Encoding(txt) <- "UTF-8"
  pages <- strsplit(paste(txt, collapse = "\n"), "\f", fixed = TRUE)[[1]]
  pages <- stringi::stri_trans_nfkc(pages)
  pages <- stringi::stri_trans_char(pages, "\uA78F\u2013\u2014\u2212\u2010\u2011", "\u00B7-----")
  lapply(pages, function(p) {
    l <- if (nzchar(p)) strsplit(p, "\r?\n")[[1]] else character()
    trimws(gsub("\\s+", " ", l, perl = TRUE))
  })
}

#' page and line for every line but the page numbers. blank lines kept as ""
pdf_lines <- function(pdf) {
  pages <- pdf_pages(pdf, layout = TRUE)
  doc <- data.frame(page = rep(seq_along(pages), lengths(pages)), line = unlist(pages))
  doc[!grepl("^[0-9]+$", doc$line, perl = TRUE), ]
}

nums <- function(line) {
  w <- split_words(line)
  w <- w[grepl("[0-9]", w)]
  sub("[.,;:)\\]}]+$", "", gsub("^[(\\[{]+|[(\\[{]+$", "", w, perl = TRUE), perl = TRUE)
}

#' a line mostly of numbers. a sentence like the end of a caption is not one
is_row <- function(line) {
  n <- length(nums(line))
  n >= 2 && n >= 0.4 * length(split_words(line)) && !endsWith(line, ".")
}

#' first run of rows in region through gaps of up to 2 other lines
#' a row repeated on a continuation page is kept once
block <- function(region) {
  rows <- character()
  gap  <- 0
  for (l in region[nzchar(region)]) {
    if (is_row(l)) {
      rows <- union(rows, l)
      gap  <- 0
    } else if (length(rows)) {
      gap <- gap + 1
      if (gap > 2) break
    }
  }
  as.character(unlist(lapply(rows, nums)))
}

#' page and row numbers in order for each table and figure
#' table rows follow the caption and figure labels come before it
displays <- function(doc) {
  i    <- seq_len(nrow(doc))
  k    <- grep("^(Appendix )?(Table|Figure) [0-9]+[:.]", doc$line, perl = TRUE)
  name <- sub("^((Appendix )?(Table|Figure) [0-9]+)[:.].*", "\\1", doc$line[k], perl = TRUE)
  out  <- list()
  for (j in seq_along(k)) {
    page <- doc$page[k[j]]
    if (grepl("Table", name[j], fixed = TRUE)) {
      end    <- if (j < length(k)) k[j + 1] - 1 else nrow(doc)
      region <- doc$line[i > k[j] & i <= min(end, k[j] + 119)]
    } else {
      start  <- if (j > 1) k[j - 1] + 1 else 1
      region <- rev(doc$line[i >= start & i < k[j] & doc$page == page])
    }
    out[[name[j]]] <- list(page = page, numbers = block(region))
  }
  out
}

#' x minus y as multisets. counts in order of first appearance in x
minus <- function(x, y) {
  u <- unique(x)
  n <- as.vector(table(factor(x, u))) - as.vector(table(factor(y[y %in% u], u)))
  setNames(n[n > 0], u[n > 0])
}

#' csv w/ quotes only where needed as the stored files have them
write_rows <- function(d, file) {
  d[] <- lapply(d, \(v) if (is.logical(v)) ifelse(v, "True", "False") else as.character(v))
  q <- function(v) ifelse(grepl("[\",\r\n]", v), paste0("\"", gsub("\"", "\"\"", v), "\""), v)
  body <- if (nrow(d)) do.call(paste, c(lapply(d, q), sep = ",")) else character()
  con <- file(file, "wb")
  writeLines(enc2utf8(c(paste(names(d), collapse = ","), body)), con, sep = "\r\n", useBytes = TRUE)
  close(con)
}

a  <- pdf_lines(stored)
b  <- pdf_lines(render)
da <- displays(a)
db <- displays(b)
six <- function(x, from) paste(head(x[seq_along(x) >= from], 6), collapse = " ")
compare <- do.call(rbind, lapply(names(da), function(name) {
  na <- da[[name]]$numbers
  nb <- if (is.null(db[[name]])) character() else db[[name]]$numbers
  m  <- min(length(na), length(nb))
  from <- c(which(na[seq_len(m)] != nb[seq_len(m)]), m + 1)[1]
  same <- identical(na, nb)
  data.frame(display = name, page = da[[name]]$page,
             render_page = if (is.null(db[[name]])) "" else as.character(db[[name]]$page),
             numbers = length(na), render_numbers = length(nb), identical = same,
             first_difference = if (same) "" else paste(six(na, from), "|", six(nb, from)))
}))
write_rows(compare, file.path(out, "render_compare.csv"))

# every number in the text w/ its pages
na_all <- lapply(a$line, nums)
nb_all <- lapply(b$line, nums)
xa <- as.character(unlist(na_all))
xb <- as.character(unlist(nb_all))
pa <- rep(a$page, lengths(na_all))
pb <- rep(b$page, lengths(nb_all))
only_a <- minus(xa, xb)
only_b <- minus(xb, xa)
pages_of <- function(v, x, pg) vapply(v, \(s) paste(sort(unique(pg[x == s])), collapse = " "), "")
diff <- data.frame(number = c(names(only_a), names(only_b)),
                   `in` = rep(c("stored PDF", "render"), c(length(only_a), length(only_b))),
                   times = unname(c(only_a, only_b)),
                   pages = unname(c(pages_of(names(only_a), xa, pa), pages_of(names(only_b), xb, pb))),
                   check.names = FALSE)
write_rows(diff, file.path(out, "render_numbers.csv"))

#' words of each page in reading order. words split at a line end are joined and hyphens
#' between letters dropped so a word reads the same wherever the line breaks
pdf_words <- function(pdf) lapply(pdf_pages(pdf, layout = FALSE), function(page) {
  joined <- character()
  for (l in page[nzchar(page) & !grepl("^[0-9]+$", page, perl = TRUE)]) {
    n <- length(joined)
    if (n && grepl("[A-Za-z]-$", joined[n], perl = TRUE) && grepl("^[a-z]", l, perl = TRUE)) {
      joined[n] <- paste0(joined[n], l)
    } else {
      joined <- c(joined, l)
    }
  }
  gsub("(?<=[A-Za-z])-(?=[A-Za-z])", "", split_words(paste(joined, collapse = " ")), perl = TRUE)
})

letters_of <- function(w) sort(as.character(unlist(strsplit(w, ""))), method = "radix")
same_letters <- function(x, y) identical(letters_of(x), letters_of(y))

# reordered = same words in other order (eg figure labels). URL split = URL or DOI broken elsewhere
# page break elsewhere = same letters once the next or previous page is added
wa <- pdf_words(stored)
wb <- pdf_words(render)
text <- do.call(rbind, lapply(seq_len(max(length(wa), length(wb))), function(p) {
  x <- if (p <= length(wa)) wa[[p]] else character()
  y <- if (p <= length(wb)) wb[[p]] else character()
  near <- intersect(c(p - 1, p + 1), seq_len(min(length(wa), length(wb))))
  verdict <- if (identical(x, y)) "identical" else
    if (identical(sort(x, method = "radix"), sort(y, method = "radix"))) "reordered" else
    if (same_letters(x, y)) "URL split elsewhere" else
    if (any(vapply(near, \(q) same_letters(c(x, wa[[q]]), c(y, wb[[q]])), TRUE))) "page break elsewhere" else
    "differs"
  sa <- minus(x, y)
  sb <- minus(y, x)
  data.frame(page = p, words = length(x), render_words = length(y), verdict = verdict,
             only_stored = paste(rep(names(sa), sa), collapse = " "),
             only_render = paste(rep(names(sb), sb), collapse = " "))
}))
write_rows(text, file.path(out, "render_text.csv"))

cat(sprintf("displays: %d, numbers identical in %d\n", nrow(compare), sum(compare$identical)))
for (r in which(!compare$identical))
  cat(sprintf("  %s (p. %s): %s\n", compare$display[r], compare$page[r], compare$first_difference[r]))
cat(sprintf("numbers in the document: stored PDF %d, render %d. %d only in the stored PDF, %d only in the render\n",
            length(xa), length(xb), sum(only_a), sum(only_b)))
v <- unique(text$verdict)
cat(sprintf("pages: stored PDF %d, render %d. %s\n", length(wa), length(wb),
            paste(table(factor(text$verdict, v)), v, collapse = ", ")))
