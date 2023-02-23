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

deltas <- list("LBX074LA" = 5,
               "LBX099LA" = 5,
               "LBX118LA" = 5,
               "LBX138LA" = 5,
               "LBX153LA" = 5,
               "LBX170LA" = 5,
               "LBX180LA" = 5,
               "LBX187LA" = 5,
               "LBX194LA" = 5,
               "LBXD03LA" = 5,
               "LBXD05LA" = 5,
               "LBXD07LA" = 5,
               "LBXF03LA" = 5,
               "LBXF04LA" = 5,
               "LBXF05LA" = 5,
               "LBXF08LA" = 5,
               "LBXHXCLA" = 5,
               "LBXPCBLA" = 5)

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
  n_folds = 20,
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
