library("hesim")
library("data.table")
library("flexsurv")
library("sets")
library("plyr")
library("MASS")
library("randomForest")
library("datasets")
library("dplyr")
library("tidyr")
library("synthpop")
library("survival")

source("../src/datasets_setup.r")
source("../src/generate_syn.r")
source("../src/distance_evaluation.r")

mean_corr = function(a, b){
  # mean pairwise correlation between two correlation matrices a and b
  return (mean(abs((a - b)[lower.tri(a - b)])))
}

fit_cox = function(data, clock="reset"){
  # fitting transition models for each transition for a semi-markov model(reset)
  # or inhomogeneous markov model (forward)
  # returns coefficients of a cox regression model
  
  n_trans = 4
  
  labs = c("Healthy-> Sick" = 1, 
           "Healthy -> Death" = 2,
           "Sick -> Healthy" = 3,
           "Sick -> Death" = 4)
  
  cox_fits = vector(length = n_trans, mode = "list")
  if (clock=='reset'){
    # adjust age so that it is age at the last transition, not age at entry
    
    data$age = data$age + data$Tstart
    cox_fits = vector(length = n_trans, mode = "list")
    for (i in 1:length(cox_fits)){
      cox_fits[[i]] = coxph(Surv(years, status) ~ treatment_grp 
                                  + age +female,
                                  data = data, 
                                  subset = (trans == i),
                            x=TRUE) 
    }
  }
  if (clock=='forward'){
    for (i in 1:length(cox_fits)){
      cox_fits[[i]] = coxph(Surv(Tstart, Tstop, status) ~ treatment_grp 
                                  + age +female,
                                  data = data, 
                                  subset = (trans == i)) 
    }
  }
  return (cox_fits)
}

iou_score = function(ia, ib){
  # calculates the interval overlap utility score
  # the intervals are vectors on the form c(lower, upper)
  iou = ((min(ia[2], ib[2])-max(ia[1], ib[1]))/(2*(ib[2]-ib[1])) 
         + (min(ia[2], ib[2])-max(ia[1], ib[1]))/(2*(ia[2]-ia[1])))
  return(iou)
}
