library(stringr)
library(zoo)
library(imputeTS)
library(here)

nhanes_data <- readRDS(here("sandbox/NHANES/output/NHANES_data.RDS"))

nhanes_data_measles <- nhanes_data[!is.na(nhanes_data$measles),]
nhanes_data_measles_pfas <- nhanes_data_measles[!is.na(nhanes_data_measles$pfas_pfdoa_ng_ml),]

nhanes_data_measles_na_thresh <- nhanes_data_measles[ , colSums(is.na(nhanes_data_measles)) < nrow(nhanes_data_measles) ]

nhanes_data_measles_na_thresh <- as.data.frame(unclass(nhanes_data_measles_na_thresh),stringsAsFactors=TRUE)

pfas <- c("pfas_pfhxs_ng_ml", "pfas_me_pfosa_acoh_ng_ml", "pfas_pfdea_ng_ml", "pfas_pfna_ng_ml", "pfas_pfua_ng_ml", "pfas_pfdoa_ng_ml")
pfas_deltas <- list("pfas_pfhxs_ng_ml" = 1, "pfas_me_pfosa_acoh_ng_ml" = 1, "pfas_pfdea_ng_ml" = 1, "pfas_pfna_ng_ml" = 1, "pfas_pfua_ng_ml" = 1, "pfas_pfdoa_ng_ml" = 1)

nhanes_data_measles_na_thresh <- nhanes_data_measles_na_thresh[complete.cases(nhanes_data_measles_na_thresh[ , c("pfas_pfhxs_ng_ml")]), ]


outcome <- "measles"

covariates <- c("age_screen_years", "gender", "race", "educ_level",
                "marital_status", "alc_gm_1999", "cotinine_ng_ml", "bmi_kg_m2",
                "fam_poverty_ratio", "systolic_1", "diastolic_1", "birth_country",
                "avg_daily_physical_act", "muscle_training",
                "vigorous_intense_30_days", "fasting_glucose_mg_dl",
                "two_year_exam_weight", "two_year_interview_weight")

w <- nhanes_data_measles_na_thresh[, covariates]
a <- nhanes_data_measles_na_thresh[, pfas]
y <- nhanes_data_measles_na_thresh$measles


                          nhanes_results <- SuperNOVA(
                            w = w,
                            a = a,
                            y = y,
                            delta = pfas_deltas,
                            n_folds = 3,
                            num_cores = 6,
                            family = "continuous",
                            quantile_thresh = 0,
                            seed = 294580
                          )

nhanes_results_neg$`Pooled TMLE Mixture Results` %>%
  dplyr::filter(Proportion_Folds >= 0.7)

mixture_plots <- plot_mixture_results(
  v_intxn_results = nhanes_results_neg$`V-Specific Mix Results`,
  hjust = 0.8)
mixture_plots$`cadmium-thallium`


qcomp <- qgcomp(Y~X1*X2*X3*X4*X5*X6*X7+Z+Z2+Z3, expnms=c(paste("X", seq(1,7), sep = "")),
              data = niehs_data,q=4, degree = 2, B=10)
plot(qcomp)


