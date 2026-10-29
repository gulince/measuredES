# run_simulations.R
# Runs both simulation studies in sequence and saves the results.
#
# Run from the package root (setwd() there first). Takes ~3-4 hours on
# 8 cores with R = 5000; for a quick test, set R_REPS <- 100 in
# 01_simulation_main.R (~5 minutes).

for (pkg in c("metafor")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("measuredES", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::load_all(".")
} else {
  library(measuredES)
}

if (!dir.exists("analysis/results")) {
  dir.create("analysis/results", recursive = TRUE)
}

message("\n=== Study 1 & 2: point estimation and CI performance ===")
message("Started: ", format(Sys.time(), "%H:%M:%S"))
source("analysis/01_simulation_main.R")
message("Done: ", format(Sys.time(), "%H:%M:%S"))

message("\n=== Study 3: meta-analytic simulation ===")
message("Started: ", format(Sys.time(), "%H:%M:%S"))
source("analysis/02_simulation_meta.R")
message("Done: ", format(Sys.time(), "%H:%M:%S"))

message("\nAll simulations finished; results are in analysis/results/.")
