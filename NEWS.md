# measuredES NEWS

## measuredES 0.1.0 (2026-07-03)

Initial release.

### Features

* `mes_d()` and `mes_g()`: Cohen's d and Hedges' g corrected for
  measurement unreliability. Both accept raw vectors, summary
  statistics (`means`, `sds`, `ns`), or a `y ~ group` formula.
* Three confidence-interval methods: the Hunter-Schmidt analytical
  approximation, Bonett's (2008) heteroscedasticity-robust variance
  with a t critical value, and the BCa bootstrap with a percentile
  fallback.
* `disattenuate()`: the bare Spearman (1904) correction, for d-family
  effect sizes and correlations.
* `mes_ci()`: recompute the confidence interval of an existing `mes`
  object with a different method; the confidence level is inherited
  from the object rather than silently reset.
* `mes_meta()`: convert one or more `mes` objects into a
  `metafor`-ready data frame (`yi`, `vi`, `sei`, `ni`).
* `mes_report()`: APA-7 inline text (English and Turkish) and tabular
  output.
* S3 methods for the `mes` class: `print()`, `summary()`, `plot()`
  (ggplot2 when available, base graphics otherwise), and
  `as.data.frame()`.

### Scope notes

* The package deliberately covers the d-family (and, via
  `disattenuate()`, correlations) only. Corrections for
  variance-explained indices (eta-squared, omega-squared) and IRT-based
  reliability were considered during development and dropped: the
  methodological literature is too thin to support them as defaults.
* All confidence-interval methods treat the reliability coefficient as
  a known constant; when it is itself estimated, the intervals are
  narrower than the full uncertainty warrants (see the vignette).
* `disattenuate(type = "r", ryy = NULL)` warns that the one-sided
  correction treats the second variable as error-free.
