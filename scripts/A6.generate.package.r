##########################################################################################
## 
##  Utils to generate an R package and the documentation associated
##  R package creation
##  https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html
##
##########################################################################################

# Once you’ve got your documentation completed, you can simply run:
devtools::document()

## Create object x in the data/ folder
usethis::use_data(x, overwrite = T)

## All packages listed in Depends, Imports and LinkingTo must be installed, 
## and their version requirements must be met, otherwise your package can’t be checked. 
## An easy way to install any missing or outdated dependencies is to run 
devtools::install_deps(dependencies = TRUE)

## We can use the devtools::load_all() function which will load your R package into memory
## exposing all the functions and data that we highlighted above. However as soon as you close your R session, 
## the package will no longer be available.
devtools::load_all() 

## To actually install your package, you use the devtools::install() function which installs your R package 
## into your R system library. Then you will be able to load up your package with library("myfirstpackage")
## Along with all the data that comes with the package!
devtools::install() 

# check automatically builds and checks a source package, using all known best practices. 
# check_built checks an already built package.
devtools::check(manual=TRUE)

## Create package manual in pdf, with all the functions
devtools::build_manual()

# You can build all vignettes from the console with devtools::build_vignettes(), but this is rarely useful. 
# Instead use devtools::build() to create a package bundle with the vignettes included. 
# RStudio’s “Build & reload” does not build vignettes to save time. 
# Similarly, devtools::install_github() (and friends) will not build vignettes by default 
# because they’re time consuming and may require additional packages. 
#You can force building with devtools::install_github(build_vignettes = TRUE). 
#This will also install all suggested packages.
devtools::build_vignettes() 
library(SBW)
browseVignettes("SBW")  # no funciona
