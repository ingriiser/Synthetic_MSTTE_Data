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

real_data = function(censored = TRUE, status = TRUE) {
  # censored=TRUE keep all censored patients
  # status = TRUE keep the status= 0 rows
  data = mstate3_exdata$transitions
  
  # sort by time
  data = data[order(data$patient_id, data$Tstart),]
  
  # make Tstart and Tstop consistent (hesim adds 1/365.25 to Tstop compared to mstate)
  # + 0.000000001 to prevent negative survival times
  data$Tstop = data$Tstop - 1 / 365.25 + 0.000000001
  data$years = data$Tstop - data$Tstart
  
  # rename strategy_id
  col_names = colnames(data)
  col_names[9] = 'treatment_grp'
  colnames(data) = col_names
  
  if (isFALSE(status)) {
    data = data[!data$status == 0,]
  }
  
  if (isFALSE(censored)) {
    # remove censored patients
    set_patient_id = unlist(as.set(data$patient_id))
    
    # patients that are still alive at the time cut-off
    censored_patients = c()
    
    for (id in set_patient_id) {
      sub = data[data$patient_id == id & data$status == 1,]
      nrows = nrow(sub)
      # if only censored transitions
      if (nrows == 0) {
        censored_patients = c(censored_patients, id)
      }
      else{
        # if the last transition is not a death
        if (sub$to[nrows] != 3) {
          censored_patients = c(censored_patients, id)
        }
      }
    }
    data = data[!(data$patient_id %in% censored_patients),]
  }
  # change treatment_grp to 0 and 1
  data["treatment_grp"][data["treatment_grp"] == 1] = 0
  data["treatment_grp"][data["treatment_grp"] == 2] = 1
  
  set_patient_id = unlist(as.set(data$patient_id))
  
  starting_state = c()
  for (id in set_patient_id) {
    sub = data[data$patient_id == id,]
    nrows = nrow(sub)
    state_1 = sub[1, 'from']
    starting_state = c(starting_state, rep(state_1, nrows))
  }
  # add column with starting state to data
  data = cbind(data, starting_state)
  
  return (data)
}

dat_patients = function(data) {
  # returns a data.frame where each row holds the constant information about
  # each patient
  
  patients_dup = data[, c('patient_id',
                          'treatment_grp',
                          'age',
                          'female',
                          'starting_state')]
  patients = patients_dup[!duplicated(patients_dup),]
  
  return(patients)
}

train_test_set = function(data, p = 0.2) {
  # divides the data into a test and training set
  # p is a value between 0 and 1 and denotes the proportion of the data
  # used in the test set
  set_patient_id = unlist(as.set(data$patient_id))
  n_patients = length(set_patient_id)
  
  n_test = floor(n_patients * p)
  n_train = n_patients - n_test
  
  test_patients = sample(set_patient_id, n_test)
  training_patients = setdiff(set_patient_id, test_patients)
  
  test_set = data[is.element(data$patient_id, test_patients),]
  training_set = data[is.element(data$patient_id, training_patients),]
  
  return(list(training_set, test_set))
}

make_trajectories = function(dataset) {
  # converting the data to vectors containing state trajectories
  # only works for uncensored observations
  # returns a list with a vector of the trajectories for each patient in dataset
  
  set_patient_id = unlist(as.set(dataset$patient_id))
  n_patients = length(set_patient_id)
  
  #convert to days
  dataset$Tstop = round(dataset$Tstop * 365.25)
  dataset$Tstart = round(dataset$Tstart * 365.25)
  
  trajectories = vector('list', n_patients)
  for (i in 1:n_patients) {
    # only include uncensored transitions
    subset = dataset[dataset$patient_id == set_patient_id[i] &
                       dataset$status == 1,]
    max_days = max(subset$Tstop) + 1
    # +1 to include day 0
    t = rep(3, max_days + 1)
    # day 0
    t[0] = subset[1, 'Tstart']
    for (j in 1:nrow(subset)) {
      # +1 to include day0
      # Day Tstop is the first day of next state
      t[subset[j, 'Tstart'] + 1:subset[j, 'Tstop']] = subset[j, 'from']
      trajectories[[i]] = t
    }
  }
  return(trajectories)
}