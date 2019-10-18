library(rio)
library(plyr)
library(tidyverse)
#devtools::install_github("FredHasselman/casnet")
library(casnet)

library(devtools)
library(usethis)

usethis::use_description()
usethis::use_roxygen_md()
usethis::use_pipe()
usethis::use_tidy_eval()
usethis::use_tidy_description()
usethis::use_tidy_style()
#usethis::use_github()
#usethis::use_tidy_github()
usethis::use_testthat()

usethis::use_dev_package("FredHasselman/casnet")
usethis::use_package("plyr")
usethis::use_package("dplyr")
usethis::use_package("DescTools")


usethis::use_package_doc()
usethis::use_gpl3_license(name = c("Fred Hasselman; Wouter Oomens"))
usethis::use_pkgdown()


tools::Rdindex(RdFiles = "~/Documents/GitHub/randseqR/man",outFile = "INDEX")
devtools::document(roclets=c('rd', 'collate', 'namespace','vignette'))
devtools::build_vignettes()
devtools::install(build_vignettes = TRUE)
pkgdown::build_site(lazy = FALSE, preview = TRUE)

