library(stringr)
library(zoo)
library(imputeTS)
library(here)
library(SuperNOVA)

nhanes_data <- readRDS(here("sandbox/NHANES/output/NHANES_data.RDS"))

nhanes_data_measles <- nhanes_data[!is.na(nhanes_data$measles), ]
nhanes_data_measles_pfas <- nhanes_data_measles[!is.na(nhanes_data_measles$pfas_pfdoa_ng_ml), ]

nhanes_data_measles_na_thresh <- nhanes_data_measles[, colSums(is.na(nhanes_data_measles)) < nrow(nhanes_data_measles)]

nhanes_data_measles_na_thresh <- as.data.frame(unclass(nhanes_data_measles_na_thresh), stringsAsFactors = TRUE)

pfas <- c("pfas_pfhxs_ng_ml", "pfas_me_pfosa_acoh_ng_ml", "pfas_pfdea_ng_ml") #"pfas_pfna_ng_ml", "pfas_pfua_ng_ml", "pfas_pfdoa_ng_ml")
pfas_deltas <- list("pfas_pfhxs_ng_ml" = 1, "pfas_me_pfosa_acoh_ng_ml" = 1, "pfas_pfdea_ng_ml" = 1) #"pfas_pfna_ng_ml" = 1, "pfas_pfua_ng_ml" = 1, "pfas_pfdoa_ng_ml" = 1)

nhanes_data_measles_na_thresh <- nhanes_data_measles_na_thresh[complete.cases(nhanes_data_measles_na_thresh[, c("pfas_pfhxs_ng_ml")]), ]


outcome <- "measles"

covariates <- c(
  "age_screen_years", "gender", "race", "educ_level",
   "cotinine_ng_ml", "bmi_kg_m2",
  "fam_poverty_ratio",
  "avg_daily_physical_act", "muscle_training",
  "vigorous_intense_30_days", "fasting_glucose_mg_dl",
  "two_year_exam_weight", "two_year_interview_weight"
)

w <- nhanes_data_measles_na_thresh[, covariates]
a <- nhanes_data_measles_na_thresh[, pfas]
y <- nhanes_data_measles_na_thresh$measles


nhanes_results <- SuperNOVA(
  w = w,
  a = a,
  y = y,
  delta = pfas_deltas,
  n_folds = 2,
  num_cores = 6,
  family = "continuous",
  quantile_thresh = 0,
  seed = 294580
)

saveRDS(
  object = nhanes_results,
  file = here("sandbox/NHANES/output/", paste0("SuperNOVA_", "sim", ".rds"))
)
