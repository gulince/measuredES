# render_manuscript.R
# Concatenates the numbered manuscript files (00-07) and renders them to
# a single Word document, manuscript/manuscript.docx.
#
# Run from the package root:
#   source("manuscript/render_manuscript.R")
#
# Requires pandoc. RStudio bundles it; outside RStudio, point the
# RSTUDIO_PANDOC environment variable at a pandoc installation first.

files <- list.files("manuscript",
                    pattern = "^0[0-8]_.*[.]Rmd$",
                    full.names = TRUE)
files <- sort(files)

combined <- file.path("manuscript", "_manuscript_combined.Rmd")
parts <- vapply(files, function(f) {
  paste(readLines(f, encoding = "UTF-8"), collapse = "\n")
}, character(1))
writeLines(paste(parts, collapse = "\n\n"), combined, useBytes = FALSE)

rmarkdown::render(combined,
                  output_file = "manuscript.docx",
                  output_dir  = "manuscript",
                  quiet       = TRUE)

unlink(combined)
message("Wrote manuscript/manuscript.docx")
