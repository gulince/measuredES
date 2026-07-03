# Manuscript

R Markdown sources for the measuredES manuscript. Target outlet:
*Educational and Psychological Measurement*, Methodological Studies
category.

Files compile in numeric order (00-07) into a single Word document.
Tables and figures are read from `analysis/results/` and
`analysis/figures/`; if the package's statistical code changes, rerun
`analysis/run_simulations.R` and
`analysis/03_analysis_tables_figures.R` before recompiling.

Before submission:

- Add an APA-7 CSL file (e.g., `apa.csl` from the Zotero style
  repository) and the journal's Word template to the YAML header of
  `00_titlepage_abstract.Rmd`.
- APA-7 style check.
- EPM reproducibility checklist and open data/code statement (the
  GitHub repository).
- Word count check (target about 9,000 words).
