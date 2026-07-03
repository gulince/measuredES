# Manuscript

R Markdown sources for the measuredES manuscript.

The numbered files (00-08) compile in order into a single Word
document. Two build scripts are provided; run either from the package
root:

- `render_manuscript.R` produces `manuscript.docx` (full version,
  with author information).
- `render_blinded.R` produces `manuscript_blinded.docx` (author
  identity and repository URLs masked for double-anonymized review)
  and `title_page.docx` (the separate title page).

Tables and figures are read from `analysis/results/` and
`analysis/figures/`; if the package's statistical code changes, rerun
`analysis/run_simulations.R` and
`analysis/03_analysis_tables_figures.R` before recompiling.
Citations use `../paper.bib` with the APA-7 CSL style (`apa.csl`).

Before submission:

- Apply the journal's Word template to the YAML header of
  `00_titlepage_abstract.Rmd` if one is required.
- Complete the journal's reproducibility checklist and open
  data/code statement (the GitHub repository).
- Verify the word count stated on the title page.
