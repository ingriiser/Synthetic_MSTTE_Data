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

generate_full = function(training_set,
                         k,
                         clock = "reset",
                         censored = FALSE) {
  # whole process of generating synthetic data
  # k: number of synthetic patients
  # clock: 'reset' or 'forward'
  # censored: if TRUE, then the censored transitions are included in the returned
  # dataframe of synthetic data
  
  # transition matrix
  tmat = rbind(c(NA, 1, 2),
               c(3, NA, 4),
               c(NA, NA, NA))
  colnames(tmat) = rownames(tmat) = c("Healthy", "Sick", "Death")
  
  training_patients = dat_patients(training_set)
  patients_syn = gen_patients(training_patients, k = k)
  wei_fits = fit_weibull(training_set, clock)
  
  transition_times_states_syn = gen_sequence_man(patients_syn, wei_fits,
                                                 tmat, clock = clock)
  synthetic_data = seq_to_dataset(transition_times_states_syn,
                                  training_set,
                                  patients_syn,
                                  censored = censored)
  return (synthetic_data)
}

gen_patients = function(patients, k) {
  # wrapper function for synthpop that generates k synthetic patients based
  # on a patient training data set, which is returned from dat_patients
  
  visit = c('age', 'female', 'treatment_grp', 'starting_state')
  syn_data = patients[, visit]
  synthetic = syn(
    data = syn_data,
    method = c('sample', 'cart', 'cart', 'cart'),
    k = k,
    visit.sequence = visit,
    minnumlevels = 3,
    print.flag = FALSE
  )
  
  patients_synthetic = synthetic$syn
  
  # adding a patient_id
  patients_syn = cbind(c(1:k), patients_synthetic)
  colnames(patients_syn)[1] = 'patient_id'
  
  return (patients_syn)
}

fit_weibull = function(data, clock = "reset") {
  # modified from hesim example code
  # fitting transition models for each transition for a semi-markov model(reset)
  # or inhomogeneous markov model (forward)
  # returns coefficients of WPHR models
  
  n_trans = 4
  
  labs = c(
    "Healthy-> Sick" = 1,
    "Healthy -> Death" = 2,
    "Sick -> Healthy" = 3,
    "Sick -> Death" = 4
  )
  
  wei_fits = vector(length = n_trans, mode = "list")
  if (clock == 'reset') {
    # adjust age so that it is age at the last transition, not age at entry
    
    data$age = data$age + data$Tstart
    wei_fits = vector(length = n_trans, mode = "list")
    for (i in 1:length(wei_fits)) {
      wei_fits[[i]] = flexsurvreg(
        Surv(years, status) ~ treatment_grp
        + age + female,
        data = data,
        subset = (trans == i) ,
        dist = "weibullPH"
      )
    }
  }
  if (clock == 'forward') {
    for (i in 1:length(wei_fits)) {
      wei_fits[[i]] = flexsurvreg(
        Surv(Tstart, Tstop, status) ~ treatment_grp
        + age + female,
        data = data,
        subset = (trans == i) ,
        dist = "weibullPH"
      )
    }
  }
  wei_fits = flexsurvreg_list(wei_fits)
  return (wei_fits)
  
}

gen_sequence_man = function(patients_syn,
                            wei_fits,
                            tmat = NULL,
                            clock = 'reset') {
  # manual synthesising of patient sequences based on the hesim procedure
  # returns a list with length = nrow(patients_syn)
  # each element contains a list consisting of two vectors of the same length
  # the first vector is the transition times of each transition from t=0 to
  # time of death
  # the second vector gives the transitions made at the times given in the first vector
  # first element is 0 to denote that no transition is made at t=0
  # use the clock that was used fitting wei_fits
  
  if (is.null(tmat)) {
    tmat = rbind(c(NA, 1, 2),
                 c(3, NA, 4),
                 c(NA, NA, NA))
    colnames(tmat) = rownames(tmat) = c("Healthy", "Sick", "Death")
  }
  
  # the number of transitions in tmat is 4
  coeff = vector('list', 4)
  for (c in  1:4) {
    coeff[[c]] = wei_fits[[c]]$coefficients
  }
  # 2 is number of transient states
  valid_trans = vector('list', 2)
  for (d in 1:2) {
    valid_trans[[d]] = tmat[d, !is.na(tmat[d, ])]
  }
  #next_ind = FALSE
  times_trans = vector('list', nrow(patients_syn))
  for (n in 1:nrow(patients_syn)) {
    r = patients_syn[n, "starting_state"]
    covariates = c(1, unlist(patients_syn
                             [n, c("treatment_grp", "age", "female")]))
    trans_times = c(0)
    # 'transition' into starting state is marked as 0
    state_trans = c(0)
    while (r != 3) {
      time_temp = c()
      for (tr in valid_trans[[r]]) {
        scale_temp = exp(coeff[[tr]][2:5] %*% covariates)
        shape_temp = exp(coeff[[tr]][1])
        if (clock == 'reset') {
          r_time = rweibullPH(1, shape = shape_temp, scale = scale_temp)
        }
        if (clock == 'forward') {
          # sample from a left-truncated distribution
          t_lower = tail(trans_times, n = 1)
          # find the quantile of T
          quantile_T = pweibullPH(tail(trans_times, n = 1),
                                  shape = shape_temp,
                                  scale = scale_temp)
          # random sample from the uniform dist to get a quantile
          u = runif(1, min = quantile_T, max = 1)
          r_time = qweibullPH(u, shape = shape_temp, scale = scale_temp)
        }
        time_temp = c(time_temp, r_time)
      }
      new_trans = valid_trans[[r]][which.min(time_temp)]
      # rounding UP to nearest day
      new_time = ceiling(min(time_temp) * 365.25) / 365.25
      if (clock == 'reset') {
        # since the time is reset each transition, we need to adjust for
        # to total time
        #new_time = tail(trans_times,n=1)+new_time + 1/365.25
        new_time_diff = new_time
        new_time = tail(trans_times, n = 1) + new_time_diff
      }
      if (clock == 'forward' & new_time == Inf) {
        # inf new_time is Inf, then all possible transition times are Inf
        # happes if the left truncation becomes pushed too far into the tail
        # then the quantile will be rounded up to 1, and r_time==Inf
        # then so much time has passed that the last event will be adjusted
        # to a death
        # trans 1->trans 2
        # trans 3 -> trans 4
        state_trans = c(head(state_trans,-1), tail(state_trans, 1) + 1)
        r = 3
        # next patient (break out of while)
        next
        
      }
      trans_times = c(trans_times, new_time)
      # finds the "to"-state of transition new_trans
      r = which(tmat == new_trans, arr.ind = TRUE)[2]
      state_trans = c(state_trans, new_trans)
      # update age SHOULD THIS JUST BE FOR TIME RESET??
      if (clock == 'reset') {
        # update the age covariate so that it refers to the age at the last transition
        covariates[3] = covariates[3] + new_time_diff
      }
    }
    times_trans[[n]] = list(trans_times, state_trans)
  }
  return (times_trans)
}

seq_to_dataset = function(times_trans,
                          data,
                          patients_synthetic,
                          censored = FALSE) {
  # transform times_trans into a data set on the same form as the real data
  # with columns patient_id, from, to, trans, Tstart, Tstop,
  # status, years, treatment_grp, age, female, starting_state
  # times_trans is an element returned by gen_sequence_man
  # patients_synthetic is an element return from gen_patients
  # if censored=TRUE the data set will also contain the censored transitions
  
  from_to_mat = rbind(c(1, 1, 2, 2), c(2, 3, 1, 3))
  rownames(from_to_mat) = c('From', 'To')
  tmat = rbind(c(NA, 1, 2),
               c(3, NA, 4),
               c(NA, NA, NA))
  colnames(tmat) = rownames(tmat) = c("Healthy", "Sick", "Death")
  
  # divide by 2 because of times_trans is a list of 2
  # subtract the number of patients, so that the death event does not get a row
  n_rows = length(unlist(times_trans)) / 2 - nrow(patients_synthetic)
  syn_data = as.data.frame(matrix(NA, nrow = n_rows, ncol = 12))
  colnames(syn_data) = colnames(data)
  # exclude the last time (death)
  syn_data['Tstart'] = unlist(lapply(times_trans, function(x) {
    head(x[[1]],-1)
  }))
  syn_data['Tstop'] = unlist(lapply(times_trans, function(x) {
    # exclude time0
    x[[1]][-1]
  }))
  syn_data['years'] = syn_data['Tstop'] - syn_data['Tstart']
  # exclude first 'transition'
  syn_data['trans'] = unlist(lapply(times_trans, function(x) {
    x[[2]][-1]
  }))
  
  # finding the number of rows per patient
  l_patients = lapply(times_trans, function(x)
    sapply(x[1], length))
  # 1 less element because of death
  l_patients = as.numeric(l_patients) - 1
  # fill in constant values
  syn_data['patient_id'] = rep(patients_synthetic$patient_id, l_patients)
  syn_data['age'] = rep(patients_synthetic$age, l_patients)
  syn_data['female'] = rep(patients_synthetic$female, l_patients)
  syn_data['treatment_grp'] = rep(patients_synthetic$treatment_grp, l_patients)
  syn_data['starting_state'] = rep(patients_synthetic$starting_state, l_patients)
  syn_data['status'] = rep(1, n_rows)
  
  syn_data['from'] = as.vector(unlist(sapply(unname(syn_data['trans']), function(x) {
    from_to_mat['From', x]
  })))
  syn_data['to'] = as.vector(unlist(sapply(unname(syn_data['trans']), function(x) {
    from_to_mat['To', x]
  })))
  
  if (censored) {
    #add censored data
    cen_df = syn_data
    cen_df$status = rep(0, n_rows)
    for (i in 1:n_rows) {
      d = syn_data[i, 'from']
      valid_trans = as.vector(tmat[d, !is.na(tmat[d, ])])
      # transition is the other possible transition
      new_trans = valid_trans[valid_trans != syn_data[i, 'trans']]
      cen_df[i, 'to'] = from_to_mat["To", new_trans]
      cen_df[i, 'trans'] = new_trans
    }
    # combine censored and uncensored data
    syn_data = as.data.frame(mapply(rbind, syn_data, cen_df))
  }
  return(syn_data)
}

discretize_time = function(time) {
  # time is a float
  interval = 1 / 365.25
  time_convered = ceiling(time / (interval)) * interval
  return (time_convered)
}
