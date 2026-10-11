# Simulation Study Guide

This directory contains the simulation studies reported in the
`measuredES` manuscript.

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

# Worked example (manuscript Section 3.4; requires psychmeta)
source("analysis/04_example_reanalysis.R")

# Bootstrap sub-study (manuscript Section 5.2.4)   # ~10 minutes
source("analysis/05_simulation_bootstrap.R")

# Review-response sub-studies: rxx_se coverage, unbalanced and
# heteroscedastic designs, severe skewness, two-sided correlation
# correction (manuscript Sections 5.2.4-5.2.6, 5.4)   # ~5 minutes
source("analysis/06_simulation_substudies.R")

# Meta-analytic pooling of corrected correlations via mes_meta()
# (manuscript Section 5.5)   # ~7 minutes
source("analysis/07_simulation_meta_correlation.R")
```

## Outputs

### Tables (`results/`)
- `simulation_main.csv` : Studies 1 & 2 raw results (192 rows)
- `simulation_meta.csv` : Study 3 raw results (12 rows)
- `table1_bias.csv`     : Bias summary (manuscript Table 1)
- `table2_coverage.csv` : CI coverage summaries (reported in Section 5.2)
- `table3_meta.csv`     : Meta-analytic comparison (manuscript Table 2)
- `substudy_meta_correlation.csv` : Pooled corrected-correlation coverage
  via `mes_meta()` (manuscript Table 7)

### Figures (`figures/`)
- `fig1_bias.{png,pdf}`    : Bias comparison (Figure 1)
- `fig2_coverage.{png,pdf}`: CI coverage probabilities (Figure 2)
- `fig3_meta_bias.{png,pdf}`: Meta-analytic bias (Figure 3)

## Reproducibility

Every simulation condition is seeded individually inside its own worker
process, so the results do not depend on the number of cores, the
operating system, or the order in which the parallel backend schedules
the conditions. Rerunning a script with the same R and package versions
reproduces its results exactly. The seeds are:

| Script | Seed for condition / cell `id` |
|---|---|
| `01_simulation_main.R` (Studies 1-2) | `20261000 + condition_id` |
| `02_simulation_meta.R` (Study 3) | `20262000 + condition_id` |
| `05_simulation_bootstrap.R` | `20260000 + cell_id` |
| `06_simulation_substudies.R` (A, B, C, D) | `20270000`, `20280000`, `20290000`, `20300000` `+ cell_id` |
| `07_simulation_meta_correlation.R` | `20263000 + condition_id` |
| `08_knapp_hartung_check.R` | single sequential stream, `set.seed(20261011)` |

The reported results were produced with R 4.1.2 (metafor 3.0.2). Set the
environment variable `MES_CORES` to choose the number of worker processes
in scripts 01, 02 and 07.

## Design notes

- R = 5,000 replications per condition. The Monte Carlo standard error
  of a coverage estimate near .95 is then about 0.003 (about 0.007 for
  Study 3 at R = 1,000), small enough that the reported differences
  between methods are not simulation noise.
- Bootstrap-BCa CI is too costly for the full factorial (5,000
  replications x 999 resamples across 192 conditions), so it is
  evaluated in the focused sub-study of
  `05_simulation_bootstrap.R` (16 cells, paired against the
  analytical intervals on the same replicates).
