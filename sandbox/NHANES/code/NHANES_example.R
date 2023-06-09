library(stringr)
library(zoo)
library(imputeTS)
library(here)
library(devtools)
load_all()
# File paths
input_path <- here("sandbox/NHANES/output/NHANES_data.RDS")
output_path <- here("sandbox/NHANES/output/", paste0("SuperNOVA_", "NHANES", ".rds"))

# Load data
nhanes_data <- readRDS(input_path)

# Variables
metals <- c("barium", "cadmium", "cobalt", "cesium", "molybdenum", "lead", "antimony", "thallium","tungsten")
metal_deltas <- setNames(rep(1, length(metals)), metals)
outcome <- "asthma"

covariates <- c(
  "age_screen_years",
  "gender",
  "race",
  "educ_level",
  "cotinine_ng_ml",
  "bmi_kg_m2",
  "alc_gm_1999",
  "fam_poverty_ratio",
  "avg_daily_physical_act",
  "muscle_training",
  "vigorous_intense_30_days",
  "two_year_exam_weight",
  "two_year_interview_weight",
  "birth_country",
  "caff_mg_1999"
)

mediators <- c(
  "mean_telomere",
  "wbc_count",
  "sd_telomere",
  "monocyte_perc",
  "neutrophils_perc",
  "c_reactive_p",
  "alpha_carotene",
  "beta_carotene",
  "vitamin_a",
  "vitamin_e",
  "trans_lycopene",
  "lutein_and_zeaxanthin"
)

# Preprocessing data
df <- nhanes_data[, c(outcome, covariates, metals, mediators), drop = FALSE]
df <- df[complete.cases(df[, c(outcome, "mean_telomere", metals)]), ]

# Retain only columns where less than 20% data is missing
df <- df[, colSums(is.na(df)) < .2 * nrow(df)]

# Imputation function
impute_mean_or_mode <- function(x) {
  if (is.numeric(x)) {
    ifelse(is.na(x), mean(x, na.rm = TRUE), x)
  } else {
    ifelse(is.na(x), names(which.max(table(x))), x)
  }
}

# Impute missing values
df_imputed <- data.frame(lapply(df, impute_mean_or_mode))

# Discretize exposures
discretize <- function(x) {
  if (length(unique(x)) < 10) {
    return(x)  # or other action
  }
  jittered_x <- jitter(x, factor = 1e-5)  # Increase the jitter factor
  cut(jittered_x,
      breaks = quantile(jittered_x, probs = seq(0, 1, 1/10), na.rm = TRUE),
      labels = FALSE, include.lowest = TRUE)
}

df_imputed[, metals] <- lapply(df_imputed[, metals], discretize)

# Run SuperNOVA
nhanes_results <- SuperNOVA(
  w = df_imputed[, covariates],
  a = df_imputed[, metals],
  z = df_imputed[, mediators],
  y = ifelse(df_imputed[, outcome] ==1, 1, 0),
  deltas = metal_deltas,
  n_folds = 10,
  num_cores = 20,
  outcome_type = "binary",
  mediator_type = "continuous",
  quantile_thresh = 0,
  seed = 294580,
  exposure_quantized = TRUE,
  mediator_quantized = FALSE,
  var_sets = NULL
)

# Save results
saveRDS(nhanes_results, output_path)
