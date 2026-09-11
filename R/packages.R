## R/packages.R --------------------------------------------------------------
## The only file in R/ with side effects: attaches every package the shared
## functions rely on and pins the dplyr verbs. Source it first, before any
## other R/ file, from every run script and qmd:
##   source(here::here("R", "packages.R"))

library(balnet)      # balnet(), cv.balnet(), balweights(): estimators_ipw.R, cv run scripts
library(glmnet)      # cv.glmnet(): estimators_ipw.R, wz run scripts
library(sbw)         # sbw(): estimators_cv.R (also attaches Matrix, quadprog, slam)
library(dplyr)       # summarise.R, plots.R (attached after sbw so dplyr::summarise wins)
library(tidyr)       # summarise.R, plots.R
library(ggplot2)     # plots.R
library(patchwork)   # (l | r) panels, plots.R
library(gt)          # tables, plots.R

# pin dplyr verbs against masking (MASS, stats, ...)
select <- dplyr::select
filter <- dplyr::filter
