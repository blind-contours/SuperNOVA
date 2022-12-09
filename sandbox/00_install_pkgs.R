# be sure to set this env variable by `export R_LIBDIR=/path/to/your/R/libs`
r_libdir <- Sys.getenv("R_LIBDIR")

# set user-specific package library
if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths(r_libdir)
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# set CRAN mirror
options(repos = structure(c(CRAN = "https://cran.rstudio.com/")))

# lie to pkgbuild, as per Jeremy
pkgbuild:::cache_set("has_compiler", TRUE)

# from CRAN
install.packages(c("here", "tidyverse", "remotes", "future", "future.apply",
                   "doFuture", "foreach", "doRNG", "data.table", "devtools",
                   "Rsolnp", "nnls", "glmnet", "Rcpp", "origami", "hal9001",
                   "speedglm"),
                 lib = r_libdir)

# use remotes to install from GitHub
remotes::install_github(c("tlverse/sl3@master",
                          "tlverse/tmle3@master"),
                        lib = r_libdir)

# update all packages
update.packages(ask = FALSE, lib.loc = r_libdir)
