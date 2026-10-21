# path floor of the Wong & Chan runs in one place. wc_helpers.R passes it to every balnet fit
# floor = last lambda reached (lam_end per rep). a path that cant reach 1e-4 ends early

#' path floor = largest |SMD| left at the end of every path (balnet max.imbalance)
wc_floor <- 1e-4

#' solver settings of every fit in a rep. the pilot fit for a cell's grid keeps balnet defaults
wc_solver <- list(maxit = 1e4, tol = 1e-5)
