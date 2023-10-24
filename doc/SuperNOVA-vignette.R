## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----figure, echo=FALSE, out.width='100%', fig.align='center'-----------------
library(knitr)
include_graphics("Biometrics_Flow_Chart.png")

## ----deltas, message=FALSE, warning=FALSE-------------------------------------
deltas <- c("M1" = 1, "M2" = 2.3, "M3" = 1.4)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(data.table)
library(dplyr)
library(kableExtra)
library(SuperNOVA)

seed <- 325911

## ----NIEHS example------------------------------------------------------------
data("NIEHS_data_1", package = "SuperNOVA")

## ----NIEH Nodes---------------------------------------------------------------
NIEHS_data_1$W <- rnorm(nrow(NIEHS_data_1), mean = 0, sd = 0.1)
w <- NIEHS_data_1[, c("W", "Z")]
a <- NIEHS_data_1[, c("X1", "X2", "X3", "X4", "X5", "X6", "X7")]
y <- NIEHS_data_1$Y

## ----run SuperNOVA NIEHS data, eval = TRUE------------------------------------
deltas <- list(
  "X1" = 1, "X2" = 1, "X3" = 1,
  "X4" = 1, "X5" = 1, "X6" = 1, "X7" = 1
)

ptm <- proc.time()

NIEH_results <- SuperNOVA(
  w = w,
  a = a,
  y = y,
  deltas = deltas,
  estimator = "tmle",
  fluctuation = "standard",
  n_folds = 2,
  outcome_type = "continuous",
  quantile_thresh = 0,
  verbose = TRUE,
  parallel = FALSE,
  parallel_type = "sequential",
  num_cores = 2,
  seed = seed,
  adaptive_delta = TRUE
)

proc.time() - ptm

indiv_shift_results <- NIEH_results$`Indiv Shift Results`
em_results <- NIEH_results$`Effect Mod Results`
joint_shift_results <- NIEH_results$`Joint Shift Results`

## ----individual results-------------------------------------------------------
indiv_shift_results$X7 %>%
  kableExtra::kbl(caption = "Effect Modification Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## ----joint results------------------------------------------------------------
em_results$X7Z %>%
  kableExtra::kbl(caption = "Effect Modification Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----interaction results------------------------------------------------------
joint_shift_results$X2X7 %>%
  kableExtra::kbl(caption = "Interaction Results") %>%
  kableExtra::kable_classic(full_width = T, html_font = "Cambria")

