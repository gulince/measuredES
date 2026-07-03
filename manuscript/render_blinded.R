# render_blinded.R
# Builds the two documents masked (double-anonymized) review requires:
#   manuscript/manuscript_blinded.docx  (no author identity, URLs masked)
#   manuscript/title_page.docx          (title, authors, affiliations)
#
# Run from the package root. Requires pandoc (see render_manuscript.R).

files <- sort(list.files("manuscript",
                         pattern = "^0[0-8]_.*[.]Rmd$",
                         full.names = TRUE))
parts <- vapply(files, function(f) {
  paste(readLines(f, encoding = "UTF-8"), collapse = "\n")
}, character(1))
txt <- paste(parts, collapse = "\n\n")

# ---- blinded main document ----
blind <- txt
blind <- gsub("https://github\\.com/gulince/measuredES",
              "[repository URL blinded for review]", blind)
blind <- gsub("(?s)author:.*?^date:", "date:", blind, perl = TRUE)
blind <- gsub("Fazilet Gül İnce Aracı", "[blinded for review]",
              blind, fixed = TRUE)
blind <- gsub("gulincegazi@gmail.com", "[blinded]", blind, fixed = TRUE)

combined <- file.path("manuscript", "_manuscript_combined.Rmd")
writeLines(blind, combined, useBytes = FALSE)
rmarkdown::render(combined,
                  output_file = "manuscript_blinded.docx",
                  output_dir  = "manuscript",
                  quiet       = TRUE)
unlink(combined)

# ---- separate title page ----
title_page <- '---
title: "Measurement-Corrected Effect Sizes for the Behavioral Sciences: A Simulation Study and the measuredES R Package"
author:
  - "Fazilet Gül İnce Aracı, Department of Measurement and Evaluation, Karamanoğlu Mehmet Bey University (gulincegazi@gmail.com, ORCID 0000-0001-5620-6911)"
date: "`r format(Sys.Date(), \'%B %Y\')`"
output:
  word_document: default
---

**Corresponding author:** Fazilet Gül İnce Aracı, Department of
Measurement and Evaluation, Karamanoğlu Mehmet Bey University.
E-mail: gulincegazi@gmail.com. ORCID: 0000-0001-5620-6911.

**Funding:** The authors declare no external funding for this work.

**Conflicts of interest:** None declared.

**Data availability:** All simulation code, results, and the
manuscript source are openly available in the package repository
(https://github.com/gulince/measuredES).

*Word count (excluding abstract, references, tables, and figures):*
approximately 5,500
'
tp <- file.path("manuscript", "_title_page.Rmd")
writeLines(title_page, tp, useBytes = FALSE)
rmarkdown::render(tp,
                  output_file = "title_page.docx",
                  output_dir  = "manuscript",
                  quiet       = TRUE)
unlink(tp)

message("Wrote manuscript/manuscript_blinded.docx and manuscript/title_page.docx")
