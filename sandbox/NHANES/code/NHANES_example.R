library(stringr)
library(zoo)
library(imputeTS)
library(here)
library(SuperNOVA)

nhanes_data <- readRDS(here("sandbox/NHANES/output/NHANES_data.RDS"))

nhanes_data_asthma <- nhanes_data[!is.na(nhanes_data$asthma), ]
nhanes_data_asthma_telomere <- nhanes_data_asthma[!is.na(nhanes_data_asthma$mean_telomere), ]
nhanes_data_asthma_telomere_cesium <- nhanes_data_asthma_telomere[!is.na(nhanes_data_asthma_telomere$cesium), ]

nhanes_data_asthma_telomere_cesium <- nhanes_data_asthma_telomere_cesium[, colSums(is.na(nhanes_data_asthma_telomere_cesium)) < .2 * nrow(nhanes_data_asthma_telomere_cesium)]

nhanes_data_asthma_telomere_cesium <- as.data.frame(unclass(nhanes_data_asthma_telomere_cesium), stringsAsFactors = TRUE)

metals <- c("barium", "cadmium", "cobalt", "cesium", "molybdenum", "lead", "antimony", "thallium","tungsten" )
metal_deltas <- list("barium" = 1, "cadmium" = 1, "cobalt" = 1, "cesium" = 1, "molybdenum" = 1, "lead" = 1, "antimony" = 1, "thallium" = 1, "tungsten" = 1 )

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



w <- nhanes_data_asthma_telomere_cesium[, covariates]
a <- nhanes_data_asthma_telomere_cesium[, metals]
z <- nhanes_data_asthma_telomere_cesium[, mediators]
y <- nhanes_data_asthma_telomere_cesium$asthma
y <- ifelse(y ==1, 1, 0)

a_imputed <- data.frame(apply(a, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))
z_imputed <- data.frame(apply(z, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))


impute_mean_or_mode <- function(w) {

  for (i in seq_along(w)) {
    if(is.numeric(w[[i]])) {
      w[[i]][is.na(w[[i]])] <- mean(w[[i]], na.rm = TRUE)
    } else if(is.factor(w[[i]])) {
      mode_val <- names(which.max(table(w[[i]])))
      w[[i]][is.na(w[[i]])] <- mode_val
    }
  }

  return(w)
}

w_imputed <- impute_mean_or_mode(w)

# Define a function to jitter and discretize a variable into quantiles
discretize <- function(x, n_quantiles) {
  x <- jitter(x, factor = 1e-10)  # Jitter the data to break ties
  cut(x, breaks = quantile(x, probs = seq(0, 1, 1/n_quantiles)),
      labels = FALSE, include.lowest = TRUE)
}

# Apply this function to each column in the data frame
a_imputed_discretized <- as.data.frame(lapply(a_imputed, discretize, n_quantiles = 10))  # Replace 4 with the desired number of quantiles

nhanes_results <- SuperNOVA(
  w = w_imputed,
  a = a_imputed_discretized,
  z = z_imputed,
  y = y,
  deltas = metal_deltas,
  n_folds = 10,
  num_cores = 20,
  outcome_type = "binary",
  quantile_thresh = 0,
  seed = 294580,
  exposure_quantized = TRUE,
  mediator_quantized = FALSE,
  var_sets = NULL
)

saveRDS(
  object = nhanes_results,
  file = here("sandbox/NHANES/output/", paste0("SuperNOVA_", "NHANES", ".rds"))
)
