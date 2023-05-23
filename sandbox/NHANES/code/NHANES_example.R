library(stringr)
library(zoo)
library(imputeTS)
library(here)
library(SuperNOVA)

nhanes_data <- readRDS(here("sandbox/NHANES/output/NHANES_data.RDS"))

nhanes_data_asthma <- nhanes_data[!is.na(nhanes_data$asthma), ]
nhanes_data_asthma_telomere <- nhanes_data_asthma[!is.na(nhanes_data_asthma$mean_telomere), ]
nhanes_data_asthma_telomere_cesium <- nhanes_data_asthma_telomere[!is.na(nhanes_data_asthma_telomere$cesium), ]

nhanes_data_asthma_telomere_cesium <- nhanes_data_asthma_telomere_cesium[, colSums(is.na(nhanes_data_asthma_telomere_cesium)) < nrow(nhanes_data_asthma_telomere_cesium)]

nhanes_data_asthma_telomere_cesium <- as.data.frame(unclass(nhanes_data_asthma_telomere_cesium), stringsAsFactors = TRUE)

metals <- c("barium", "cadmium", "cobalt", "cesium", "molybdenum", "lead", "antimony", "thallium","tungsten" )
metal_deltas <- list("barium" = 1, "cadmium" = 1, "cobalt" = 1, "cesium" = 1, "molybdenum" = 1, "lead" = 1, "antimony" = 1, "thallium" = 1, "tungsten" = 1 )

outcome <- "asthma"

covariates <- c(
  "age_screen_years", "gender", "race", "educ_level",
   "cotinine_ng_ml", "bmi_kg_m2",
  "fam_poverty_ratio",
  "avg_daily_physical_act", "muscle_training",
  "vigorous_intense_30_days", "fasting_glucose_mg_dl",
  "two_year_exam_weight", "two_year_interview_weight"
)

w <- nhanes_data_asthma_telomere_cesium[, covariates]
a <- nhanes_data_asthma_telomere_cesium[, metals]
z <- nhanes_data_asthma_telomere_cesium[, "mean_telomere"]
y <- nhanes_data_asthma_telomere_cesium$asthma
y <- ifelse(y ==1, 1, 0)

nhanes_results <- SuperNOVA(
  w = w,
  a = a,
  z = z,
  y = y,
  delta = metal_deltas,
  n_folds = 2,
  num_cores = 6,
  outcome_type = "binary",
  quantile_thresh = 0,
  seed = 294580
)

saveRDS(
  object = nhanes_results,
  file = here("sandbox/NHANES/output/", paste0("SuperNOVA_", "sim", ".rds"))
)
