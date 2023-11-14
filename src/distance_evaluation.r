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

source("../src/datasets_setup.r")
source("../src/generate_syn.r")

comp_traj = function(t1, t2) {
  # compare two trajectories using the normalised hamming distance
  # the trajectories are elements of a list returned by make_trajectories(dataset)
  # max distance = 1, min distance=0
  score = t1 == t2
  return(1 - sum(score) / length(t1))
}


find_closest = function(pat_id,
                        data1,
                        data2,
                        age_int = Inf,
                        same_sex = FALSE,
                        same_t = FALSE,
                        w = 0.5,
                        nndr = TRUE) {
  # for a patient in data1 with patient id pat_id, find the closest patient in
  # data2. The search can be restricted to patients of the same sex, treatment
  # group and/or within an age interval.
  # weight balances the distance between covariates and state trajectories
  # dist = w*dist_traj + (1-w)*cov_dist
  # if nndr==TRUE it also finds the second closest patient
  
  
  # finding the covariates of the patient
  age_p = data1[data1$patient_id == pat_id, 'age'][1]
  sex_p = data1[data1$patient_id == pat_id, 'female'][1]
  treat_p = data1[data1$patient_id == pat_id, 'treatment_grp'][1]
  
  # trajectory of the patient
  traj_p = make_trajectories(data1[data1$patient_id == pat_id, ])[[1]]
  
  if (same_sex & same_t) {
    #find the matrix consisting of all observations in data2 with the same
    # covariates as the patient above
    sub = data2[data2$treatment_grp == treat_p
                & data2$sex == sex_p
                & data2$age >= age_p - age_int / 2
                & data2$age <= age_p + age_int / 2, ]
  }
  if (same_sex & isFALSE(same_t)) {
    sub = data2[data2$sex == sex_p
                & data2$age >= age_p - age_int / 2
                & data2$age <= age_p + age_int / 2, ]
  }
  if (same_t & isFALSE(same_sex)) {
    sub = data2[data2$treatment_grp == treat_p
                & data2$age >= age_p - age_int / 2
                & data2$age <= age_p + age_int / 2, ]
  }
  else{
    sub = data2[data2$age >= age_p - age_int / 2
                & data2$age <= age_p + age_int / 2, ]
  }
  
  if (nrow(sub) == 0) {
    age_int = age_int * 1.5
    if (age_int > 100) {
      msg = 'There are no elements in data2 that are close to the patient with the given patient_id.
Consider also measuring distances between patients of different sexes or different treatment groups.'
      stop(msg)
    }
    print(paste0(
      'No match. Trying again by increasing the age interval to ',
      age_int,
      '.'
    ))
    
    find_closest(pat_id, data1, data2, age_int = age_int, w = w)
  }
  
  # ids of the patients we will measure the distance to
  set_ids_data2 = unlist(as.set(sub$patient_id))
  
  # find the distance to each
  distances = sapply(set_ids_data2, function(x) {
    # covariate distance
    age_x = sub[sub$patient_id == x, 'age'][1]
    sex_x = sub[sub$patient_id == x, 'female'][1]
    treat_x = sub[sub$patient_id == x, 'treatment_grp'][1]
    
    dist_cov = norm_cov_dist(age_p, sex_p, treat_p,
                             age_x, sex_x, treat_x, data1)
    
    
    # trajectory distance
    traj_x = make_trajectories(data2[data2$patient_id == x, ])[[1]]
    
    if (length(traj_p) != length(traj_x)) {
      # add extra 3s to the shortest vector
      max_l = max(length(traj_p), length(traj_x))
      traj_p = c(traj_p, rep(3, max_l - length(traj_p)))
      traj_x = c(traj_x, rep(3, max_l - length(traj_x)))
    }
    dist_traj = comp_traj(traj_p, traj_x)
    
    #total distance
    
    dist = w * dist_cov + (1 - w) * dist_traj
    return(dist)
  })
  
  if (nndr) {
    ind_sort = order(distances)
    ind_min = ind_sort[1]
    # distance to closest patient
    min_dist = distances[ind_min]
    ind_2nd_min = ind_sort[2]
    # distance to second closest patient
    min_dist_2nd = distances[ind_2nd_min]
    if (ind_min == 0) {
      nndr = 0
    }
    else{
      nndr = min_dist / min_dist_2nd
    }
    
    return (list(
      min_dist,
      set_ids_data2[ind_min],
      min_dist_2nd,
      set_ids_data2[ind_2nd_min],
      nndr
    ))
  }
  
  else{
    ind_min = which.min(distances)
    min_dist = distances[ind_min]
    return(list(min_dist, set_ids_data2[ind_min]))
  }
}

norm_cov_dist = function(age1, sex1, treat1, age2, sex2, treat2, data) {
  # normalize the covariate distance measures according to a data set
  # returns a normalized total distance measure
  
  # mean sd normalization to have most values within [-1,1]
  age_mean = mean(data$age)
  age_sd = sd(data$age)
  sex_mean = mean(data$female)
  sex_sd = sd(data$female)
  treat_mean = mean(data$treatment_grp)
  treat_sd = sd(data$treatment_grp)
  
  age_dist = ((age1 - age_mean) / age_sd - (age2 - age_mean) / age_sd) ^
    2
  sex_dist = ((sex1 - sex_mean) / sex_sd - (sex2 - sex_mean) / sex_sd) ^
    2
  treat_dist = ((treat1 - treat_mean) / treat_sd - (treat2 - treat_mean) /
                  treat_sd) ^ 2
  
  dist = sqrt(age_dist + sex_dist + treat_dist) / sqrt(3 * 4)
  return(dist)
}

wilcox = function(a, b, alternative = "two.sided") {
  # wrapper function for the exact wilcox test in the package coin.
  # a and b are two vectors
  # returns the p.value
  data_df <- data.frame(Value = c(a, b),
                        Group = factor(c(
                          rep("Group1", length(a)), rep("Group2", length(b))
                        )))
  res = wilcox_test(
    Value ~ Group,
    data = data_df,
    distribution = "exact",
    alternative = alternative
  )
  p = pvalue(res)
  return (p)
  
  
}

evaluate_closeness = function(data1,
                              data2,
                              age_int = Inf,
                              same_sex = FALSE,
                              same_t = FALSE,
                              weight = 0.5,
                              nndr = TRUE) {
  # for each element in data1 find the distance to the closest element in data2
  # return a vector with the distances
  
  set_ids = unlist(as.set(data1$patient_id))
  if (identical(data1, data2)) {
    # intra distance
    if (nndr) {
      result = lapply(set_ids, function(x) {
        # removing the current patient id from data2
        data2_new = data2[data2$patient_id != x, ]
        return(
          find_closest(
            x,
            data1,
            data2_new,
            age_int = age_int,
            same_sex = same_sex,
            same_t = same_t,
            w = weight,
            nndr = TRUE
          )[c(1, 5)]
        )
      })
      return (result)
    }
    
    else{
      result = lapply(set_ids, function(x) {
        # removing the current patient id from data2
        data2_new = data[data$patient_id != x, ]
        return(
          find_closest(
            x,
            data1,
            data2_new,
            age_int = age_int,
            same_sex = same_sex,
            same_t = same_t,
            w = weight,
            nndr = FALSE
          )[[1]]
        )
      })
      
      return (result)
    }
  }
  else{
    # inter distance
    if (nndr) {
      result = lapply(set_ids, function(x) {
        find_closest(
          x,
          data1,
          data2,
          age_int = age_int,
          same_sex = same_sex,
          same_t = same_t,
          w = weight,
          nndr = TRUE
        )[c(1, 5)]
      })
      return (result)
    }
    
    else{
      result = lapply(set_ids, function(x) {
        find_closest(
          x,
          data1,
          data2,
          age_int = age_int,
          same_sex = same_sex,
          same_t = same_t,
          w = weight,
          nndr = FALSE
        )[[1]]
      })
      
      return (result)
    }
  }
  
}