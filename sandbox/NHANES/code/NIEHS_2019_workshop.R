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

deltas <- list("LBX074LA" = 2,
               "LBX099LA" = 2,
               "LBX118LA" = 2,
               "LBX138LA" = 2,
               "LBX153LA" = 2,
               "LBX170LA" = 2,
               "LBX180LA" = 2,
               "LBX187LA" = 2,
               "LBX194LA" = 2,
               "LBXD03LA" = 2,
               "LBXD05LA" = 2,
               "LBXD07LA" = 2,
               "LBXF03LA" = 2,
               "LBXF04LA" = 2,
               "LBXF05LA" = 2,
               "LBXF08LA" = 2,
               "LBXHXCLA" = 2,
               "LBXPCBLA" = 2)

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
  adaptive_delta = FALSE
)

saveRDS(
  object = nhanes_results,
  file = here("sandbox/NHANES/output", paste0("SuperNOVA_", "nhanes_prime", ".rds"))
)
