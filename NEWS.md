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
* Optional `rxx_se` argument in `mes_d()` and `mes_g()`: supplying
  the standard error of the reliability estimate adds the
  first-order delta-method term to the analytical variances, so the
  interval reflects that reliability was estimated rather than
  known. `mes_ci()` and `mes_meta()` honor it.
* `disattenuate()`: the bare Spearman (1904) correction, for d-family
  effect sizes and correlations.
* `mes_r()`: correlation corrected for one- or two-sided
  unreliability, with a Fisher-z confidence interval whose endpoints
  are divided by the attenuation factor (exact under known
  reliabilities; endpoints truncated to [-1, 1]).
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
* By default the confidence intervals treat the reliability
  coefficient as a known constant; when it is itself estimated, pass
  `rxx_se` to propagate that uncertainty (the bootstrap always holds
  `rxx` fixed).
* `disattenuate(type = "r", ryy = NULL)` warns that the one-sided
  correction treats the second variable as error-free.
