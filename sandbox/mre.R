library(here)
library(sl3)
library(tmle3)

source(here("01_setup_data.R"))
devtools::load_all(here())

data <- sim_data(400)

cv_folds <- 5

## setup learners

# continuous outcome learners
fglm_contin_lrnr <- Lrnr_glm_fast$new()
lasso_contin_lrnr <- Lrnr_glmnet$new(alpha = 1, family = "gaussian",
                                     nfolds = cv_folds)

# binary outcome learners
fglm_binary_lrnr <- Lrnr_glm_fast$new(family = binomial())
lasso_binary_lrnr <- Lrnr_glmnet$new(alpha = 1, family = "binomial",
                                     nfolds = cv_folds)

sl_contin <- Lrnr_sl$new(
  learner = list(fglm_contin_lrnr, lasso_contin_lrnr),
  metalearner = Lrnr_nnls$new()
)

sl_binary <- Lrnr_sl$new(
  learner = list(fglm_binary_lrnr, lasso_binary_lrnr),
  metalearner = Lrnr_nnls$new()
)

node_list <- list(
  W = c("W1", "W2", "W3"),
  A = "A",
  Z = c("Z1", "Z2", "Z3"),
  Y = "Y"
)

learner_list <- list(
  Y = sl_contin,
  A = sl_binary
)

tmle_spec_NDE_error <- tmle_NDE(
  e_learners = sl_binary,
  psi_Z_learners = sl_contin
)

tmle3(tmle_spec_NDE_error, data, node_list, learner_list)

## Error in column_names[names(new_col_map)] <- new_col_map : 
##            'names' attribute [11] must be the same length as the vector [8]
## Failed on chain
## Error in self$compute_step() : 
##            Error in column_names[names(new_col_map)] <- new_col_map : 
##                       'names' attribute [11] must be the same length as the vector [8]

sessionInfo()

## R version 3.6.3 (2020-02-29)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.4 LTS

## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     

## other attached packages:
##  [1] tmle3mediate_0.0.1 testthat_2.3.2     forcats_0.4.0      stringr_1.4.0     
##  [5] dplyr_0.8.5        purrr_0.3.3        readr_1.3.1        tidyr_1.0.2       
##  [9] tibble_3.0.1       ggplot2_3.3.0      tidyverse_1.3.0    tmle3_0.1.7       
## [13] sl3_1.3.7          here_0.1          

## loaded via a namespace (and not attached):
##  [1] nlme_3.1-142         fs_1.3.1             usethis_1.6.1       
##  [4] lubridate_1.7.4      devtools_2.3.0       progress_1.2.2      
##  [7] httr_1.4.1           rprojroot_1.3-2      tools_3.6.3         
## [10] backports_1.1.6      R6_2.4.1             DBI_1.0.0           
## [13] colorspace_1.4-1     withr_2.2.0          tidyselect_1.0.0    
## [16] prettyunits_1.1.1    processx_3.4.2       compiler_3.6.3      
## [19] glmnet_3.0-2         cli_2.0.2            rvest_0.3.5         
## [22] xml2_1.2.2           desc_1.2.0           scales_1.1.0        
## [25] checkmate_2.0.0      nnls_1.4             mvtnorm_1.1-0       
## [28] randomForest_4.6-14  callr_3.4.3          speedglm_0.3-2      
## [31] digest_0.6.25        pkgconfig_2.0.3      htmltools_0.4.0     
## [34] sessioninfo_1.1.1    dbplyr_1.4.2         htmlwidgets_1.5.1   
## [37] rlang_0.4.6          readxl_1.3.1         rstudioapi_0.11     
## [40] BBmisc_1.11          shape_1.4.4          visNetwork_2.0.9    
## [43] generics_0.0.2       jsonlite_1.6.1       magrittr_1.5        
## [46] delayed_0.3.0        Matrix_1.2-17        Rcpp_1.0.4.6        
## [49] munsell_0.5.0        fansi_0.4.1          abind_1.4-5         
## [52] lifecycle_0.2.0      stringi_1.4.6        MASS_7.3-51.4       
## [55] pkgbuild_1.0.7       grid_3.6.3           parallel_3.6.3      
## [58] listenv_0.8.0        crayon_1.3.4         lattice_0.20-38     
## [61] haven_2.2.0          hms_0.5.3            ps_1.3.2            
## [64] pillar_1.4.3         igraph_1.2.5         uuid_0.1-4          
## [67] origami_1.0.3        future.apply_1.5.0   codetools_0.2-16    
## [70] pkgload_1.0.2        reprex_0.3.0         glue_1.4.0          
## [73] data.table_1.12.8    remotes_2.1.1        modelr_0.1.5        
## [76] vctrs_0.2.4          foreach_1.5.0        cellranger_1.1.0    
## [79] gtable_0.3.0         rstackdeque_1.1.1    future_1.17.0       
## [82] assertthat_0.2.1     broom_0.5.2          iterators_1.0.12    
## [85] memoise_1.1.0        globals_0.12.5       ellipsis_0.3.0      
## [88] imputeMissings_0.0.3
