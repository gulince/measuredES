# measuredES: Measurement-Corrected Effect Sizes

> Effect sizes corrected for measurement error, for the behavioral sciences.

[![R-CMD-check](https://github.com/gulince/measuredES/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gulince/measuredES/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

`measuredES` computes Cohen's *d* and Hedges' *g* corrected for
measurement unreliability. If you have a reliability estimate for your
outcome measure, you can use it to get a better approximation of the
underlying effect, with less distortion from measurement noise.

Packages like `effectsize`, `MBESS`, and `psych` each cover a piece of
this, and `psychmeta` corrects *d* values for unreliability inside its
meta-analysis framework. `measuredES` handles the single-study side, in
one place: the correction itself (from raw data, summary statistics, or
a formula), three confidence interval methods (Hunter-Schmidt, Bonett,
Bootstrap BCa), APA-7 ready output, and direct `metafor` integration.

## What it does

Under classical test theory every observed score is *X = T + E* (true
score plus noise), and that noise pulls standardized effect sizes
toward zero:

```
d_observed = d_true × √(r_xx)
```

So at *r_xx* = 0.70 the observed effect is about 16% smaller than its
construct-level counterpart. `measuredES` reverses that shrinkage,
drawing on the psychometric tradition (Hunter & Schmidt, 2004) and
modern interval estimation (Bonett, 2008).

## Installation

```r
# Install from GitHub
remotes::install_github("gulince/measuredES")
```

## Quick start

```r
library(measuredES)

set.seed(2026)
control   <- rnorm(60, mean = 100, sd = 15)
treatment <- rnorm(60, mean = 108, sd = 15)

# From raw vectors
mes_d(x1 = control, x2 = treatment, rxx = 0.85)

# From summary statistics (e.g., values reported in a published article)
mes_d(means = c(100, 108), sds = c(15, 15), ns = c(60, 60), rxx = 0.85)

# Formula interface: outcome ~ group, with a data frame
df <- data.frame(
  score = c(control, treatment),
  group = rep(c("control", "treatment"), each = 60)
)
mes_d(formula = score ~ group, data = df, rxx = 0.85)
```

### A note on sign

`mes_d()` computes `d = (mean(x1) - mean(x2)) / Sp`, so if `x2` is the
group with the larger mean (a typical treatment effect), `d` comes out
negative. Swap `x1` and `x2` if you want the conventional "treatment
positive" sign.

## Functions

| Function | What it does |
|---|---|
| `mes_d()` | Cohen's *d* with reliability correction |
| `mes_g()` | Hedges' *g* with reliability correction |
| `mes_r()` | Correlation corrected for unreliability, with a Fisher-z CI |
| `disattenuate()` | The bare math behind the correction (d and r) |
| `mes_ci()` | Re-run the CI on an existing `mes` object with a different method |
| `mes_meta()` | Build a `metafor`-ready data frame (`yi`, `vi`, `sei`, `ni`) from one or more `mes` objects |
| `mes_report()` | APA-7 inline text or a one-row table for the manuscript |

For a quick flat view of one result, `as.data.frame(fit)` gives back a
single-row frame with the estimate, CI bounds, method, reliability, and
n. For meta-analysis, use `mes_meta()`: it's the one that produces the
`yi`/`vi` columns `metafor::rma()` expects.

## References

- Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. *Psychological Methods*, *13*(2), 99–109. https://doi.org/10.1037/1082-989X.13.2.99
- Charles, E. P. (2005). The correction for attenuation due to measurement error: Clarifying concepts and creating confidence sets. *Psychological Methods*, *10*(2), 206–226. https://doi.org/10.1037/1082-989X.10.2.206
- Hunter, J. E., & Schmidt, F. L. (2004). *Methods of meta-analysis: Correcting error and bias in research findings* (2nd ed.). Sage.
- Spearman, C. (1904). The proof and measurement of association between two things. *American Journal of Psychology*, *15*(1), 72–101. https://doi.org/10.2307/1412159

## Contributing and support

- **Bug reports and feature requests:** please open an issue at
  <https://github.com/gulince/measuredES/issues>, ideally with a minimal
  reproducible example.
- **Contributions:** pull requests are welcome. Please make sure
  `devtools::check()` passes and add a test for any behavior change.
- **Questions:** open a GitHub issue, or contact the maintainer at
  <gulincegazi@gmail.com>.

## License

GPL (>= 3)
