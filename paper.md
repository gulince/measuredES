---
title: 'measuredES: An R Package for Measurement-Corrected Effect Sizes in the Behavioral Sciences'
tags:
  - R
  - psychometrics
  - effect size
  - measurement error
  - reliability
  - meta-analysis
authors:
  - name: Fazilet Gül İnce Aracı
    orcid: 0000-0001-5620-6911
    affiliation: 1
  - name: Eren Halil Özberk
    orcid: 0000-0003-2136-3081
    affiliation: 2
affiliations:
 - name: Department of Measurement and Evaluation, Karamanoğlu Mehmet Bey University
   index: 1
 - name: National Foundation for Educational Research, Slough, United Kingdom
   index: 2
date: 3 July 2026
bibliography: paper.bib
---

# Summary

`measuredES` is an R package that computes effect sizes
(Cohen's *d* and Hedges' *g*) corrected for attenuation due to
measurement unreliability. The package implements psychometric correction traditions
[@hunter2004; @spearman1904] with modern interval estimation
techniques [@charles2005; @bonett2008]. Outputs are compatible with
`metafor` for meta-analysis and follow APA-7 reporting standards.

# Statement of need

The misreporting of effect sizes due to measurement error is a
documented yet under-corrected problem across the behavioral
sciences [@schmidt1996; @baugh2002]. When measures have
imperfect reliability ($r_{xx} < 1$), observed standardized effect
sizes systematically underestimate population effects:

$$d_{observed} = d_{true} \cdot \sqrt{r_{xx}}$$

This attenuation distorts power analyses, replicability assessments,
and meta-analytic syntheses. Existing R packages provide effect size
computation [`effectsize`: @bensh2020; `MBESS`: @kelley2007;
`psych`: @revelle2024] but none offer a unified framework for
measurement-corrected *d*-family effect sizes with multiple confidence
interval methods and meta-analytic output. `measuredES` fills this gap.

# Core functionality

The starting point is `disattenuate()`, which applies Spearman's
correction to any observed *d*-family value. Built on top of that are
`mes_d()` and `mes_g()`, which accept raw vectors, a formula, or the
summary statistics you pulled from a published paper.

Confidence intervals are available via three methods: the
Hunter--Schmidt analytical approximation, Bonett's (2008)
heteroscedasticity-robust variance, and the bias-corrected accelerated
bootstrap [@efron1987]. For balanced equal-variance designs the first
two agree closely; Bonett's method becomes important when group SDs
differ or the design is unbalanced.

`mes_report()` produces inline APA-7 text that can go straight into a
manuscript. `mes_meta()` converts results to the format `metafor::rma()`
expects, so moving from a single study to a meta-analysis is one function
call.

# Example

```r
library(measuredES)
set.seed(2026)
control <- rnorm(60, mean = 100, sd = 15)
treat   <- rnorm(60, mean = 108, sd = 15)
fit <- mes_d(x1 = control, x2 = treat,
             rxx = 0.78, ci_method = "hunter_schmidt")
mes_report(fit, style = "inline", language = "en")
#> corrected d = -0.61 (observed d = -0.53), 95% CI [-1.02, -0.19];
#>  reliability r_xx = 0.78
```

# Acknowledgements

The authors declare no external funding for this work.

# References
