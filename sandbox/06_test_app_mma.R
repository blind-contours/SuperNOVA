
# mma package install - probably not install for most people:

list.of.packages <- c("mma")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# packages
library(here)
library(foreach)
library(future)
library(doFuture)
library(doRNG)
library(data.table)
library(tidyverse)
library(hal9001)
library(origami)
library(sl3)
library(tmle3)
library(mma)
devtools::load_all(here())
`%notin%` <- Negate(`%in%`)

##get simulating data sources
source(here("sandbox/01_setup_data.R"))


##load data from mma for bmi study
data("weight_behavior")

##remove NA as per stochastic paper
weight_behavior_na_rm <- weight_behavior[complete.cases(weight_behavior), ] #567, 124 obs removed

##refactor exposure variable to 0,1:

weight_behavior_na_rm$sports <- ifelse(weight_behavior_na_rm$sports == 2, 1, 0)

##set up node list - same as stochastic paper
exposure <- "sports"
outcome <- "bmi"
mediators <- c("snack", "exercises", "overweigh")

##set up node list
covars <- colnames(weight_behavior_na_rm[colnames(weight_behavior_na_rm) %notin% c(exposure, outcome, mediators)])

node_list <- list(
  W = covars,
  A = exposure,
  Z = mediators,
  Y = outcome
)

#set up learners
hal_lrnr <- Lrnr_hal9001$new(max_degree = NULL,
                             n_folds = 5,
                             fit_type = "glmnet",
                             use_min = TRUE,
                             type.measure = "deviance",
                             standardize = FALSE,
                             family = "gaussian",
                             lambda.min.ratio = 1 / nrow(data),
                             nlambda = 500,
                             yolo = FALSE)


hal_contin_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5
)
hal_binary_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5,
  family = "binomial"
)

cv_hal_contin_lrnr <- Lrnr_cv$new(hal_contin_lrnr, full_fit = TRUE)
cv_hal_binary_lrnr <- Lrnr_cv$new(hal_binary_lrnr, full_fit = TRUE)


learner_list <- list(
  Y = cv_hal_contin_lrnr,
  A = cv_hal_binary_lrnr
)

##nie
tmle_spec_NIE <- tmle_NIE(
  e_learners = cv_hal_binary_lrnr,
  psi_Z_learners = cv_hal_contin_lrnr,
  max_iter = 100
)

weight_behavior_NIE <- tmle3(tmle_spec_NIE, weight_behavior_na_rm, node_list, learner_list)
weight_behavior_NIE <- weight_behavior_NIE$summary

print(weight_behavior_NIE)

##nde
tmle_spec_NDE <- tmle_NDE(
  e_learners = cv_hal_binary_lrnr,
  psi_Z_learners = cv_hal_contin_lrnr,
  max_iter = 100
)

weight_behavior_NDE <- tmle3(tmle_spec_NDE, weight_behavior_na_rm, node_list, learner_list)
weight_behavior_NDE <- weight_behavior_NDE$summary

print(weight_behavior_NDE)

###################### MMA method #####################

x = weight_behavior_na_rm[,c(covars, exposure, mediators)]

pred = weight_behavior_na_rm[, exposure]
y = data.frame(weight_behavior_na_rm[,outcome])
mediators <- weight_behavior_na_rm[, mediators]

#binary predictor
#binary y
x=weight_behavior_na_rm[,c(2:11,13:15)]
pred=weight_behavior_na_rm[,12]
y=data.frame(weight_behavior_na_rm[,1])
colnames(y)="bmi"


mma_glm<-mma(x,
             y,
             pred=pred,
             mediator=c(7,11,13),
             jointm=list(n=1,j1=c(7,11,13)),
             predref=0,
             alpha=0.4,
             alpha2=0.4,
             n2=20,
             nonlinear=FALSE)

print(mma_glm)

mma_mart <-mma(x,
              y,
              pred=pred,
              mediator=c(7,11,13),
              jointm=list(n=1,j1=c(7,11,13)),
              predref=0,
              alpha=0.4,
              alpha2=0.4,
              n =2,
              n2=1,
              nonlinear=TRUE)

print(mma_mart)

