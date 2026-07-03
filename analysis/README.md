# Simulation Study Guide

This directory contains the simulation studies for the `measuredES` manuscript
submitted to *Educational and Psychological Measurement*.

## Prerequisites

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "metafor"))
devtools::install("..")  # install measuredES from parent directory
```

## Execution order

```r
# Set working directory to the package root
setwd("path/to/measuredES")

# Studies 1 & 2: Point estimation and CI performance
source("analysis/01_simulation_main.R")   # ~3-4 hours (8-core machine)

# Study 3: Meta-analytic
source("analysis/02_simulation_meta.R")   # ~30 minutes

# Tables and figures
source("analysis/03_analysis_tables_figures.R")
```

## Outputs

### Tables (`results/`)
- `simulation_main.csv` : Studies 1 & 2 raw results (192 rows)
- `simulation_meta.csv` : Study 3 raw results (12 rows)
- `table1_bias.csv`     : Bias summary (manuscript Table 1)
- `table2_coverage.csv` : CI coverage (manuscript Table 2)
- `table3_meta.csv`     : Meta-analytic comparison (Table 3)

### Figures (`figures/`)
- `fig1_bias.{png,pdf}`    : Bias comparison (Figure 1)
- `fig2_coverage.{png,pdf}`: CI coverage probabilities (Figure 2)
- `fig3_meta_bias.{png,pdf}`: Meta-analytic bias (Figure 3)

## Reproducibility

All scripts begin with `set.seed(2026)`. Results may show minor
differences across operating systems due to parallel RNG; the number of
cores used should be reported alongside manuscript results.

## Design notes

- R = 5,000 replications per condition; sufficient for stable Monte Carlo
  estimates at EPM standards.
- Bootstrap-BCa CI is excluded from the main simulation due to
  computational cost (5,000 replications x 2,000 resamples per cell is
  prohibitive across 192 conditions).
