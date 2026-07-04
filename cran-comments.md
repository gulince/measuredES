# cran-comments

## Test environments

* Local: Windows 11, R 4.5.3
* GitHub Actions: Ubuntu (R release and R devel), Windows (R release),
  macOS (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

(The only NOTE on the built tarball comes from the CRAN incoming
feasibility check: new submission.)

## Comments

* This is a new submission (first release of the package).
* Author names in DESCRIPTION contain non-ASCII characters (Turkish);
  the package declares `Encoding: UTF-8`.
* Words the spell checker may flag in DESCRIPTION (Bonett, Hedges,
  disattenuation, psychometric, metafor, APA) are proper names and
  standard methodological terms; an `inst/WORDLIST` is included.
