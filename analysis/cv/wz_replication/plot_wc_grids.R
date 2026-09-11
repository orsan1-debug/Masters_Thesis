# plot_wc_grids.R ----------------------------------------------------------
# Writes the path grids, lambda-pick figures and lambda / RMSE tables of
# wc_grids.R to output/cv/wz_replication/ (figures/, tables/), for the noise x overlap, n x overlap
# and basis x overlap batches. The qmd sources wc_grids.R and draws the same
# figures inline; this script is for looking at them outside the render.

fig_dir <- here::here("output", "cv", "wz_replication", "figures")
source(here::here("analysis", "cv", "wz_replication", "wc_grids.R"))

# *** Noise x overlap (Sections 2.1 / 2.3 analogue) ***
sig <- c(1, 10, 30, 100); cc <- 1:4
noise <- read_grid("wc_noise_s", sig, cc)
for (o in c("A", "B")) {
  grid_fig(noise, sig, cc, "sigma", "c", o,
           title = sprintf("Wong & Chan design, ATT, model %s: noise x overlap (n = 5000, 1000 reps)", o),
           file = file.path(fig_dir, paste0("fig_wc_noise_overlap_", o, ".png")))
}
pick_fig(noise, sig, cc, "sigma", "c",
         title = "Median selected lambda vs RMSE-optimal lambda: noise axis (top), overlap axis (bottom)",
         file = file.path(fig_dir, "fig_wc_noise_overlap_lambda.png"))
write_tabs(noise, sig, cc, "sigma", "c", "wc_noise")

# *** n x overlap ***
ns <- c(1000, 2500, 5000, 10000, 20000)
on <- read_grid("wc_on_n", ns, cc)
for (o in c("A", "B")) {
  grid_fig(on, ns, cc, "n", "c", o,
           title = sprintf("Wong & Chan design, ATT, model %s: n x overlap (sigma = 1, 1000 reps)", o),
           file = file.path(fig_dir, paste0("fig_wc_n_overlap_", o, ".png")))
}
pick_fig(on, ns, cc, "n", "c",
         title = "Median selected lambda vs RMSE-optimal lambda: n axis (top), overlap axis (bottom)",
         file = file.path(fig_dir, "fig_wc_n_overlap_lambda.png"))
write_tabs(on, ns, cc, "n", "c", "wc_on")

# *** basis size x overlap ***
ks <- c(10, 20, 65, 125); cb <- c(1, 3)
basis <- read_grid("wc_basis_K", ks, cb)
for (o in c("A", "B")) {
  grid_fig(basis, ks, cb, "K", "c", o,
           title = sprintf("Wong & Chan design, ATT, model %s: basis size x overlap (n = 5000, sigma = 1, 1000 reps)", o),
           file = file.path(fig_dir, paste0("fig_wc_basis_overlap_", o, ".png")))
}
pick_fig(basis, ks, cb, "K", "c",
         title = "Median selected lambda vs RMSE-optimal lambda: basis axis (top), overlap axis (bottom)",
         file = file.path(fig_dir, "fig_wc_basis_overlap_lambda.png"))
write_tabs(basis, ks, cb, "K", "c", "wc_basis")
