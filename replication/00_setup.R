# checks packages and installs balnet 0.0.4 from pkg/ if needed. writes replication/sessionInfo.txt
#   Rscript replication/00_setup.R   the render also needs Quarto XeLaTeX and Times New Roman

# R.utils lets fread read .csv.gz. ggh4x knitr rmarkdown for the render. readr to rebuild summaries
needed <- c("here", "data.table", "R.utils", "dplyr", "tidyr", "purrr", "tibble", "glmnet",
            "ggplot2", "gt", "withr", "testthat", "Rcpp", "RcppEigen", "future", "furrr",
            "ggh4x", "knitr", "rmarkdown", "readr")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

if (!requireNamespace("balnet", quietly = TRUE) ||
    packageVersion("balnet") != "0.0.4") {
  tarball <- here::here("pkg", "balnet_0.0.4.tar.gz")
  stopifnot(file.exists(tarball))
  install.packages(tarball, repos = NULL, type = "source")   # needs a C++17 compiler
}

library(balnet); library(glmnet)
invisible(lapply(needed, loadNamespace))   # load all so sessionInfo records their versions
# stored results used R 4.5.2 (win) and 4.6.1 (mac) w/ balnet 0.0.4 and glmnet 5.0
# wc cells used glmnet 4.1.10. glmnet 4.1.8 gave the stored glmnet est to 1e-13 (smoke_test.R)
cat("R", R.version.string, "\nbalnet", as.character(packageVersion("balnet")),
    "\nglmnet", as.character(packageVersion("glmnet")), "\n")
writeLines(c(capture.output(sessionInfo()),
             paste("balnet:", packageVersion("balnet")),
             paste("glmnet:", packageVersion("glmnet"))),
           here::here("replication", "sessionInfo.txt"))
cat("written replication/sessionInfo.txt\n")
