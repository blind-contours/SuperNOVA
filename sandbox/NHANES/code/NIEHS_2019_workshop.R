library(SuperNOVA)
library(readr)
library(here)

NIEHS_2019 <- read_csv(here("sandbox/NHANES/input/NIEHS_2019.csv"))

exposures <- c("LBX074LA",
               "LBX099LA",
               "LBX118LA",
               "LBX138LA",
               "LBX153LA",
               "LBX170LA",
               "LBX180LA",
               "LBX187LA",
               "LBX194LA",
               "LBXD03LA",
               "LBXD05LA",
               "LBXD07LA",
               "LBXF03LA",
               "LBXF04LA",
               "LBXF05LA",
               "LBXF08LA",
               "LBXHXCLA",
               "LBXPCBLA")

NIEHS_2019 <- NIEHS_2019[complete.cases(NIEHS_2019[, exposures]), ]

deltas <- list("LBX074LA" = 1,
               "LBX099LA" = 1,
               "LBX118LA" = 1,
               "LBX138LA" = 1,
               "LBX153LA" = 1,
               "LBX170LA" = 1,
               "LBX180LA" = 1,
               "LBX187LA" = 1,
               "LBX194LA" = 1,
               "LBXD03LA" = 1,
               "LBXD05LA" = 1,
               "LBXD07LA" = 1,
               "LBXF03LA" = 1,
               "LBXF04LA" = 1,
               "LBXF05LA" = 1,
               "LBXF08LA" = 1,
               "LBXHXCLA" = 1,
               "LBXPCBLA" = 1)

outcome <- "TELOMEAN"

covariates <- c("LBXWBCSI",
                "LBXLYPCT",
                "LBXMOPCT",
                "LBXEOPCT",
                "LBXBAPCT",
                "LBXNEPCT",
                "male",
                "age_cent",
                "age_sq",
                "race_cat",
                "bmi_cat3",
                "ln_lbxcot",
                "edu_cat")

w <- NIEHS_2019[, covariates]
a <- NIEHS_2019[, exposures]
y <- NIEHS_2019$TELOMEAN


nhanes_results <- SuperNOVA(
  w = w,
  a = a,
  y = y,
  delta = deltas,
  n_folds = 10,
  num_cores = 20,
  family = "continuous",
  quantile_thresh = 0,
  seed = 294580,
  adaptive_delta = TRUE
)

saveRDS(
  object = nhanes_results,
  file = here("sandbox/NHANES/output", paste0("SuperNOVA_", "nhanes_prime", ".rds"))
)
