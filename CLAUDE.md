# Project rules
- Root = folder containing the .Rproj. Never use the folder name in
  paths or code; it will be renamed. Use here::here(), forward slashes,
  no setwd(), no absolute paths. Must work on Windows and Mac.
- results/ is read-only. Never delete; move to archive/ instead.
- Work in phases (inventory, move, path repair, verify, docs). One
  commit per phase. Stop and wait after each phase.
- No logic changes to estimators or DGPs without asking.
- R: tidyverse, snake_case, roxygen on every function.
- Read the r-sim-edit skill before touching R code.