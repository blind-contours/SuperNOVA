################################################################################
### Title:        NHANES Variable Processing
### Author:        Aidan
### Date Created:  4/14/2020
################################################################################



################################################################################
### I. Environment, data loading
################################################################################

`%p%` <- function(x,y) paste0(x,y)
library(tidyverse)
library(magrittr)
library(janitor)
library(readxl)
library(purrr)
library(gtools)
library(here)

### get variable list
var_list <-
  read_excel(here("sandbox/NHANES/input/Dictionary_1516.xlsx"), sheet = "selected",
             col_names = FALSE)

names(var_list) <-
  c("var_name", "var_description")


### get data dictionary
data_dic <-
  read_excel(here("sandbox/NHANES/input/Dictionary_1516.xlsx"), sheet = "total",
             col_names = TRUE)


names(data_dic) <-
  c("var_name", "var_description")


### get relevant NHANES data
rds_files <-
  list.files(here("sandbox/NHANES/input/"),pattern="data.rds")


empty_df <- c()


### manually adding some dietary variables for now.
diet_vars <- c("DR1TCAFF", "DR2TCAFF", "DRXTCAFF",
               "DR1TALCO", "DR2TACLO", "DRXTALCO",
               "DR1TPROT", "DR2TPROT", "DRXTPROT",
               "DR1TCHOL", "DR2TCHOL", "DRXTCHOL",
               "DR1TSFAT", "DR2TSFAT", "DRXTSFAT",
               "DR1TMFAT", "DR2TMFAT", "DRXTMFAT",
               "DR1TPFAT", "DR2TPFAT", "DRXTPFAT")

full_var_list <-
  append(var_list$var_name, diet_vars)


for(paths in rds_files) {

  assign(gsub("_data.rds", "", paths),

         readRDS(here("sandbox/NHANES/input/") %p% paths) %>%
           dplyr::select(matches(paste(full_var_list, collapse = "|"))) %>%
           mutate(name = gsub("_data.rds", "", paths))
  )

  if(ncol(get(gsub("_data.rds", "", paths))) <= 2) {
    empty_df <-
      append(empty_df,
             gsub("_data.rds", "", paths))}


}

### remove data if no variables of interest
rm(list=empty_df)

### some variables appear twice in the Glucose datasets.
for(g in c("B", "C", "D", "E", "F", "G", "H", "I")) {
  if("GLU_" %p% g %in% ls()) {
    assign("GLU_" %p% g,
           dplyr::select(get("GLU_" %p% g),
                  -contains("PHAFST")))
  }
}

### redundant tables
rm(COT_E)


################################################################################
### II. Data Merging
################################################################################


### Get list of data
data_list = list()
names <- list()
j=1
for(i in ls()) {

  if("data.frame" %in% class(get(i)) & !i %in% c("var_list", "data_dic")) {
    data_list[[j]] = get(i) %>% mutate(name = i)
    names[[j]] <- i
    j=j+1
  }
}

names(data_list) <- names

data_list <-
  lapply(data_list,
       function(x) x %>%
         dplyr::mutate(year_begin = dplyr::case_when(
           grepl("_A", name) ~ 1999,
           grepl("_B", name) ~ 2001,
           grepl("_C", name) ~ 2003,
           grepl("_D", name) ~ 2005,
           grepl("_E", name) ~ 2007,
           grepl("_F", name) ~ 2009,
           grepl("_G", name) ~ 2011,
           grepl("_H", name) ~ 2013,
           grepl("_I", name) ~ 2015,
           ### odd cases
           name == "LAB06"   ~ 1999,
           name == "LAB06HM"   ~ 1999,
           name == "LAB13AM" ~ 1999,
           name == "LAB10AM" ~ 1999,
           name == "DEMO"    ~ 1999,
           name == "DRXTOT"  ~ 1999,
           name == "BMX"     ~ 1999,
           name == "BPX"     ~ 1999,
           name == "PAQ"     ~ 1999,
           name == "L19_C_data"  ~ 2003,
           name == "L19_B_data"  ~ 2001,
           TRUE ~ NA_real_
         )) %>% dplyr::select(-name))


for(y in c(1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015)){
  print(y)

  year_matches_list <- data_list[unlist(lapply(data_list,
                                               function(x) sum(x$year_begin == y) > 0))]

  year_matches_list <- purrr::compact(year_matches_list)

  assign("data_" %p% y,
         purrr::reduce(year_matches_list,
                full_join,  by = c("SEQN", "year_begin")))
}


full_data <-
  bind_rows(data_1999, data_2001, data_2003, data_2005,
            data_2007, data_2009, data_2011, data_2013, data_2015)


################################################################################
### III. Variable Naming
################################################################################

full_data <-
  full_data %>%
  dplyr::select(participant_index         = SEQN,
         year_begin,
         ### survey weights
         two_year_interview_weight = WTINT2YR,
         two_year_exam_weight      = WTMEC2YR,

         ## HHREF demo vars
         hhref_age_years           = DMDHRAGE,
         hhref_gender              = DMDHRGND,
         hhref_birth_country       = DMDHRBR4,
         hhref_educ_level          = DMDHREDU,
         hhref_marital_status      = DMDHRMAR,
         hhref_spouse_educ_level   = DMDHSEDU,

         ## Study participant demo vars
         age_screen_months         = RIDAGEMN,
         age_screen_years          = RIDAGEYR,
         age_exam_months           = RIDEXAGM,
         gender                    = RIAGENDR,
         birth_country             = DMDBORN,
         birth_country_other       = DMDBORN4,
         educ_level                = DMDEDUC,
         educ_level_adults         = DMDEDUC2,
         educ_level_children       = DMDEDUC3,
         marital_status            = DMDMARTL,
         years_in_us               = DMDYRSUS,

         ## SES Variables
         fam_income                = INDFMIN2,
         fam_poverty_ratio         = INDFMPIR,
         household_income          = INDHHIN2,
         preg_status               = RIDEXPRG,
         race                      = RIDRETH1,
         race_recode               = RIDRETH3,
         ### these two seem to be redundant
         ###fam_pov_level_index       = INDFMMPC,
         ###fam_size_pov_level_ratio  = INDFMMPI,

         ## Logistical
         exam_period               = RIDEXMON,
         interview_exam_status     = RIDSTATR,
         data_release_cycle        = SDDSRVYR,
         variance_pseudo_psu       = SDMVPSU,
         variance_pseudo_stratum   = SDMVSTRA,


         ### weight variables
         interview_2year_weight    = WTINT2YR,
         exam_2year_weight         = WTMEC2YR,
         bmi_kg_m2                 = BMXBMI,
         standing_height_cm        = BMXHT,
         bmi_category_children     = BMDBMIC,
         weight_kg                 = BMXWT,


         ### Diet Chemical Variables
         energy_kcal_day1            = DR1TKCAL,
         caffeine_mg_day1            = DR1TCAFF,
         caffeine_mg_day2            = DR2TCAFF,
         caff_mg_1999                = DRXTCAFF,
         alcohol_gm_day1             = DR1TALCO,
         alcohol_gm_day2             = DR2TALCO,
         alc_gm_1999                 = DRXTALCO,
         protein_gm_day1             = DR1TPROT,
         protein_gm_day2             = DR2TPROT,
         prot_gm_1999                = DRXTPROT,
         cholesterol_mg_day1         = DR1TCHOL,
         cholesterol_mg_day2         = DR2TCHOL,
         cholesterol_mg_1999         = DRXTCHOL,
         satur_fat_gm_day1           = DR1TSFAT,
         satur_fat_gm_day2           = DR2TSFAT,
         satur_fat_gm_1999           = DRXTSFAT,
         monounsat_fat_gm_day1       = DR1TMFAT,
         monounsat_fat_gm_day2       = DR2TMFAT,
         monounsat_fat_gm_1999       = DRXTMFAT,
         polyunsat_fat_gm_day1       = DR1TPFAT,
         polyunsat_fat_gm_day2       = DR2TPFAT,
         polyunsat_fat_gm_1999       = DRXTPFAT,
         on_diet                     = DRQSDIET,
         type_of_diet                = DRQSDT1,
         fast_minutes                = PHAFSTMN,
         fast_hours                  = PHAFSTHR,
         fasting_glucose_mmol_l      = LBDGLUSI,
         fasting_glucose_mg_dl       = LBXGLU,
         triglyceride_mg_per_dl     =  LBXTR,
         triglyceride_mmol_per_l    =  LBDTRSI,
         ldl_cholesterol_mg_per_dl  =  LBDLDL,
         ldl_cholesterol_mmol_per_l =  LBDLDLSI,
         hdl_cholesterol_mg_per_dl  =  LBDHDD,
         hdl_cholesterol_mmol_l     =  LBDHDDSI,
         ### smoking proxy
         cotinine_ng_ml             =  LBXCOT,
         ### Diet Habit Variables
         diet_health = DBQ700,
         regular_milk_use = DBQ229,
         not_home_meals_past_week    = DBD895,
         fast_food_meals_past_week   = DBD900,
         ready_to_eat_foods_past_30_days = DBD905,
         frozen_meals_past_30_days = DBD910,

         ### Exercise Vars
         strength_exercise_days_past_week                = PAQ678,
         doc_told_to_exercise_12mo                       = MCQ365B,
         intense_work_days_in_week                       = PAQ610,
         mod_intensity_work_activities                   = PAQ620,
         walk_bike_10_mins_travel_days_in_week           = PAQ640,
         vigorous_intense_workout_days_in_week           = PAQ655,
         days_physically_active_more_than_hour_past_week = PAQ706,
         avg_daily_physical_act                          = PAQ180,
         vigorous_intense_30_days                        = PAD200,
         muscle_training                                 = PAD440,
         num_time_week_hard_exercise                     = PAQ560,

         ### PFAS variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ### TODO: confirm via google sheet these are PFAS of interest!
         pfas_pfhxs_ng_ml                 = LBXPFHS,
         pfas_pfhxs_comment               = LBDPFHSL,
         pfas_me_pfosa_acoh_ng_ml         = LBXMPAH,
         pfas_me_pfosa_acoh_comment       = LBDMPAHL,
         pfas_pfdea_ng_ml                 = LBXPFDE,
         pfas_pfdea_comment               = LBDPFDEL,
         pfas_pfna_ng_ml                  = LBXPFNA,
         pfas_pfna_comment                = LBDPFNAL,
         pfas_pfua_ng_ml                  = LBXPFUA,
         pfas_pfua_comment                = LBDPFUAL,
         pfas_pfdoa_ng_ml                 = LBXPFDO,
         pfas_pfdoa_comment               = LBDPFDOL,
         pfas_npfoa_ng_ml                 = LBXNFOA,
         pfas_npfoa_comment               = LBDNFOAL,
         pfas_sb_pfoa_ng_ml               = LBXBFOA,
         pfas_sb_pfoa_comment             = LBDBFOAL,
         pfas_npfos_ng_ml                 = LBXNFOS,
         pfas_npfos_comment               = LBDNFOSL,
         pfas_smpfos_ng_ml                = LBXMFOS,
         pfas_smpfos_comment              = LBDMFOSL,


         ### Diabetes variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         diabetes                         = DIQ010,

         ### Insulin variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         insulin_pmol_L                   = LBDINSI,
         insulin_uu_ml                    = LBXIN,

         ### Immune variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         wbc_count                        = LBXWBCSI,
         c_reactive_protein               = LBXCRP,
         monocyte_perc                    = LBXMOPCT,
         neutrophils_perc                 = LBXNEPCT,

         ### metal variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         barium                           = URXUBA,
         cadmium                          = URDUCD,
         cobalt                           = URXUCO,
         cesium                           = URXUCS,
         molybdenum                       = URXUMO,
         manganese                        = URXUMN,
         lead                             = URXUPB,
         antimony                         = URXUSB,
         tin                              = URXUSN,
         strontium                        = URXUSR,
         thallium                         = URXUTL,
         tungsten                         = URXUTU,

         ### telomere variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         mean_telomere                    = TELOMEAN,
         sd_telomere                      = TELOSTD,

         ### blood pressure variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         systolic_1                       = BPXSY1,
         diastolic_1                      = BPXDI1,
         systolic_2                       = BPXSY2,
         diastolic_2                      = BPXDI2,
         measles                          = LBXME,
         rubella                          = LBDRUIU,
         varicella                        = LBXVAR


         ) %>%
  dplyr::mutate(year_end = year_begin + 1)


### make sure no year mistmatch:
stopifnot(sum(full_data %$% table(data_release_cycle, year_begin) >1) ==
            length(unique(full_data$year_begin)))


full_data <-
  full_data %>%
  dplyr::select(-data_release_cycle)



### CREATE HOMAR_IR - ~~~~~~~~~~~~

full_data <- full_data %>% mutate(homa_ir = ((insulin_uu_ml)*fasting_glucose_mmol_l)/22.5)


### SUBPART feature editing and creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### TODO: Note that the variance_pseudo_stratum is broken-up by data year.


### RECODING
feature_data <-
  full_data %>%
  mutate(hhref_birth_country =
           case_when(
             hhref_birth_country == 1  ~ "USA",
             hhref_birth_country == 2  ~ "Others",
             hhref_birth_country == 77 ~ "Refused",
             hhref_birth_country == 99 ~ "Dont Know",
             TRUE ~ NA_character_),
         birth_country =
           case_when(
             birth_country == 1 ~ "USA",
             birth_country %in% c(2,3) ~ "Others",
             birth_country_other == 1 ~ "USA",
             birth_country_other == 2 ~ "Others",
             TRUE ~ NA_character_
           ),
         hhref_educ_level =
           case_when(
             hhref_educ_level == 1    ~ "Less than 9th Grade",
             hhref_educ_level == 2    ~ "9-11th Grade",
             hhref_educ_level == 3    ~ "High School Grad",
             hhref_educ_level == 4    ~ "Some College or AA Degree",
             hhref_educ_level == 5    ~ "College Graduate or Above",
             hhref_educ_level == 7    ~ "Refused",
             hhref_educ_level == 9    ~ "Dont Know",
             TRUE ~ NA_character_),
         educ_level_children =
           case_when(
             educ_level_children == 0  ~ "Never attended or Kindergarten",
             educ_level_children == 13 ~ "High School Grad",
             educ_level_children == 14 ~ "GED",
             educ_level_children == 15 ~ "More than High School",
             educ_level_children == 55 ~ "Less than 5th Grade",
             educ_level_children == 66 ~ "Less than 9th Grade",
             educ_level_children == 77 ~ "Refused",
             educ_level_children == 99 ~ "Dont Know",
             is.na(educ_level_children) ~ NA_character_,
             TRUE ~ as.character(educ_level_children)
           ),
         educ_level_adults =
           case_when(
             educ_level_adults == 1    ~ "Less than 9th Grade",
             educ_level_adults == 2    ~ "9-11th Grade",
             educ_level_adults == 3    ~ "High School Grad",
             educ_level_adults == 4    ~ "Some College or AA Degree",
             educ_level_adults == 5    ~ "College Graduate or Above",
             educ_level_adults == 7    ~ "Refused",
             educ_level_adults == 9    ~ "Dont Know",
             TRUE ~ NA_character_
           ),
         educ_level =
           case_when(
             !is.na(educ_level_adults) ~ educ_level_adults,
             TRUE ~ educ_level_children
           ),
         hhref_gender =
           case_when(
             hhref_gender == 1 ~ "Male",
             hhref_gender == 2 ~ "Female",
             ### note other code is described as "Missing"
             TRUE ~ NA_character_
           ),
         hhref_marital_status =
           case_when(
             hhref_marital_status == 1     ~ "Married",
             hhref_marital_status == 2     ~ "Widowed",
             hhref_marital_status == 3     ~ "Divored",
             hhref_marital_status == 4     ~ "Separated",
             hhref_marital_status == 5     ~ "Never married",
             hhref_marital_status == 6     ~ "Living with partner",
             hhref_marital_status == 77    ~ "Refused",
             hhref_marital_status == 99    ~ "Dont Know",
             TRUE ~ NA_character_
           ),
         hhref_spouse_educ_level =
           case_when(
             hhref_spouse_educ_level == 1    ~ "Less than 9th Grade",
             hhref_spouse_educ_level == 2    ~ "9-11th Grade",
             hhref_spouse_educ_level == 3    ~ "High School Grad",
             hhref_spouse_educ_level == 4    ~ "Some College or AA Degree",
             hhref_spouse_educ_level == 5    ~ "College Graduate or Above",
             hhref_spouse_educ_level == 7    ~ "Refused",
             hhref_spouse_educ_level == 9    ~ "Dont Know",
             TRUE ~ NA_character_
           ),
         marital_status =
           case_when(
             marital_status == 1     ~ "Married",
             marital_status == 2     ~ "Widowed",
             marital_status == 3     ~ "Divored",
             marital_status == 4     ~ "Separated",
             marital_status == 5     ~ "Never married",
             marital_status == 6     ~ "Living with partner",
             marital_status == 77    ~ "Refused",
             marital_status == 99    ~ "Dont Know",
             TRUE ~ NA_character_
           ),
         years_in_us =
           case_when(years_in_us == 1  ~ "Less than 1",
                     years_in_us == 2  ~ "Between 1 and 5",
                     years_in_us == 3  ~ "Between 5 and 10",
                     years_in_us == 4  ~ "Between 10 and 15",
                     years_in_us == 5  ~ "Between 15 and 20",
                     years_in_us == 6  ~ "Between 20 and 30",
                     years_in_us == 7  ~ "Between 30 and 40",
                     years_in_us == 8  ~ "Between 40 and 50",
                     years_in_us == 9  ~ "50 or more",
                     years_in_us == 77 ~ "Refused",
                     years_in_us == 88 | years_in_us == 99 ~ "Dont Know or Could Not Determine",
                     TRUE ~ NA_character_),
         # fam_income =
         #   case_when(fam_income == 1 ~ "0 to 5K",
         #             fam_income == 2 ~ "5K to 10K",
         #             fam_income == 3 ~ "10K to 15K",
         #             fam_income == 4 ~ "15K to 20K",
         #             fam_income == 5 ~ "20K to 25K",
         #             fam_income == 6 ~ "25K to 35K",
         #             fam_income == 7 ~ "35K to 45K",
         #             fam_income == 8 ~ "45K to 55K",
         #             fam_income == 9 ~ "55K to 65K",
         #             fam_income == 10 ~ "65K to 75K",
         #             fam_income == 12 ~ "More than 20",
         #             fam_income == 13 ~ "Under 20",
         #             fam_income == 14 ~ "75K to 100K",
         #             fam_income == 15 ~ "100K and over",
         #             fam_income == 77 ~ "Refused",
         #             fam_income == 99 ~ "Dont Know",
         #             TRUE ~ NA_character_),
         ### fam_income RECODE:
         fam_income = case_when(fam_income %in% c(1,2,3,4,13) ~ "Under 20K",
                                fam_income %in% c(5,6,7,8) ~ "20K to 55K",
                                fam_income %in% c(9,10,14) ~ "55K to 100K",
                                fam_income == 15 ~ "100K and over",
                                fam_income %in% c(12, 77, 99) ~ NA_character_,
                                TRUE ~ "Missing"),
         household_income =
           case_when(household_income == 1 ~ "0 to 5K",
                     household_income == 2 ~ "5K to 10K",
                     household_income == 3 ~ "10K to 15K",
                     household_income == 4 ~ "15K to 20K",
                     household_income == 5 ~ "20K to 25K",
                     household_income == 6 ~ "25K to 35K",
                     household_income == 7 ~ "35K to 45K",
                     household_income == 8 ~ "45K to 55K",
                     household_income == 9 ~ "55K to 65K",
                     household_income == 10 ~ "65K to 75K",
                     household_income == 12 ~ "More than 20",
                     household_income == 13 ~ "Under 20",
                     household_income == 14 ~ "75K to 100K",
                     household_income == 15 ~ "100K and over",
                     household_income == 77 ~ "Refused",
                     household_income == 99 ~ "Dont Know",
                     TRUE ~ NA_character_),
         gender = case_when(
           gender == 1 ~ "Male",
           gender == 2 ~ "Female",
           TRUE ~ NA_character_
         ),
         exam_period = case_when(
           exam_period == 1 ~ "Nov 1 to April 30",
           exam_period == 2 ~ "May 1 to Oct 31",
           TRUE ~ NA_character_
         ),
         preg_status = case_when(
           preg_status == 1 ~ "Pregnant",
           preg_status == 2 ~ "Not Pregnant",
           preg_status == 3 ~ "Cannot Ascertain",
           gender == "Male" ~ "Not Female",
           TRUE ~ NA_character_
           ),
         race = case_when(
           !is.na(race_recode) & race_recode == 6 ~ "Asian",
           !is.na(race_recode) & race_recode == 7 ~ "Other Race or Multi-Racial",
           race == 1 ~ "Mexican American",
           race == 2 ~ "Other Hispanic",
           race == 3 ~ "White",
           race == 4 ~ "Black",
           race == 5 ~ "Other Race or Multi-Racial",
           TRUE ~ NA_character_
         ),
         interview_exam_status = case_when(
           interview_exam_status == 1 ~ "Interviewed only",
           interview_exam_status == 2 ~ "Interviewed and Examined",
           TRUE ~ NA_character_
         ),
         doc_told_to_exercise_12mo = case_when(
           doc_told_to_exercise_12mo == 1 ~ "Yes",
           doc_told_to_exercise_12mo == 2 ~ "No",
           doc_told_to_exercise_12mo == 9 ~ "Dont Know"
         ),
         ### TODO: since underlying is numeric, I'll code "Don't Know" to NA where applicable
         intense_work_days_in_week = case_when(
           intense_work_days_in_week == 99 ~ NA_real_,
           TRUE ~ intense_work_days_in_week
         ),
         walk_bike_10_mins_travel_days_in_week = case_when(
           walk_bike_10_mins_travel_days_in_week == 99 ~ NA_real_,
           TRUE ~ walk_bike_10_mins_travel_days_in_week
         ),
         vigorous_intense_workout_days_in_week = case_when(
           vigorous_intense_workout_days_in_week == 99 ~ NA_real_,
           TRUE ~ vigorous_intense_workout_days_in_week
         ),
         days_physically_active_more_than_hour_past_week = case_when(
           days_physically_active_more_than_hour_past_week == 77 ~ NA_real_,
           days_physically_active_more_than_hour_past_week == 99 ~ NA_real_,
           TRUE ~ days_physically_active_more_than_hour_past_week
         ),
         mod_intensity_work_activities = case_when(
           mod_intensity_work_activities == 1 ~ "Yes",
           mod_intensity_work_activities == 2 ~ "No",
           mod_intensity_work_activities == 7 ~ "Refused",
           mod_intensity_work_activities == 9 ~ "Dont Know",
           TRUE ~ NA_character_
         ),
         ### diet variables
         on_diet = case_when(
           on_diet == 1 ~ "Yes",
           on_diet == 2 ~ "No",
           on_diet == 9 ~ "Dont Know"
         ),
         type_of_diet = case_when(
           type_of_diet == 1 ~ "Weight loss or low calorie diet",
           TRUE ~ "Other or Missing"
         ),
         bmi_category_children = case_when(
           bmi_category_children == 1 ~ "Underweight",
           bmi_category_children == 2 ~ "Normal weight",
           bmi_category_children == 3 ~ "Overweight",
           bmi_category_children == 4 ~ "Obese",
           TRUE ~ NA_character_
         ),
         ### dieting habit variables
         diet_health = case_when(
           diet_health == 1 ~ "Excellent",
           diet_health == 2 ~ "Very good",
           diet_health == 3 ~ "Good",
           diet_health == 4 ~ "Fair",
           diet_health == 5 ~ "Poor",
           diet_health == 7 ~ "Refused",
           diet_health == 9 ~ "Dont Know",
           TRUE ~ NA_character_
         ),
         regular_milk_use = case_when(
           regular_milk_use == 1 ~ "Regular Drinker Most of Life",
           regular_milk_use == 2 ~ "Never Been Regular Drinker",
           regular_milk_use == 3 ~ "Milk Drinking has Varied",
           regular_milk_use == 7 ~ "Refused",
           regular_milk_use == 9 ~ "Dont Knor",
           TRUE ~ NA_character_
         ),
         ### similarly these have numeric structure, so I NA special cases.
         not_home_meals_past_week    = case_when(
           not_home_meals_past_week > 21 ~ NA_real_,
           TRUE ~ not_home_meals_past_week
         ),
         fast_food_meals_past_week   = case_when(
           fast_food_meals_past_week > 21 ~ NA_real_,
           TRUE ~ fast_food_meals_past_week
         ),
         ready_to_eat_foods_past_30_days = case_when(
           ready_to_eat_foods_past_30_days > 150 ~ NA_real_,
           TRUE ~ ready_to_eat_foods_past_30_days
           ),
         frozen_meals_past_30_days = case_when(
           frozen_meals_past_30_days > 150 ~ NA_real_,
           TRUE ~ frozen_meals_past_30_days
         )
         ) %>%
  dplyr::select(-race_recode)


index_vars <- c("participant_index", "year_begin", "year_end")

feature_data <-
  feature_data[,c(index_vars, base::setdiff(names(feature_data),
                                         index_vars))]


### I don't include capital letters in my variables so this is an
###     effective filter
feature_data <-
  feature_data[, !grepl("[A-Z]", names(feature_data))]



### save the data.
feature_data %>%
  saveRDS(here("sandbox/NHANES/output/NHANES_data.RDS"))


rm(list = ls())


















