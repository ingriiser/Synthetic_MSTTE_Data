library("hesim")
library("data.table")
library("flexsurv")
library("sets")
library("plyr")
library("MASS")
library("randomForest")
library("datasets")
library("dplyr")
library("synthpop")
library("coin")

library("rstan")
rstan_options(auto_write = TRUE)

source("../src/datasets_setup.r")
source("../src/generate_syn.r")
source("../src/distance_evaluation.r")

fit_weibull_dp = function(data, eta, epsilon=1){
  #Works only for clock reset
  #Returns a list with the Weibull models corresponding to each transition
  # in data so that wei_fits[[c]]$coefficients returns coefficients of a
  # WPHR model
  # epsilon is the epsilon value of each transition. The total privacy budget is
  # epsilon * n_trans
  sm = stan_model("../src/dp_posterior_patient.stan")
  n_trans = max(data$trans)
  wei_fits = vector(length=n_trans, mode="list")
  for (i in 1:n_trans){
    data_i = data[data$trans==i,]
    ids_n = as.vector(table(data_i$patient_id)) # number of rows per patient
    list_input_data = list(M= length(ids_n), N=nrow(data_i), m= ids_n, 
                           X=data_i[,c(9:11)], 
                           t = data_i[,8], d = data_i[,7], 
                           eta=eta, epsilon=epsilon)
    
    wei_fits[[i]]$fit_stan = stan(
      file = "../src/dp_posterior_patient.stan", 
      data = list_input_data, 
    )
    extracted = rstan::extract(wei_fits[[i]]$fit_stan)
    ind = sample(1:length(extracted[[1]]), 1) # draw a random index
    # ONE sample from posterior
    # convert to same format as flexsurvreg
    wei_fits[[i]]$coefficients = c(shape = log(extracted$lambda[ind]), 
                                   scale = extracted$beta0[ind], 
                                   treatment_grp = extracted$beta[ind,1],
                                   age = extracted$beta[ind,2],
                                   female = extracted$beta[ind,3])
  }
  return (wei_fits)
}
