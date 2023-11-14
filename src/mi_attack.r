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

source("../src/datasets_setup.r")
source("../src/generate_syn.r")
source("../src/distance_evaluation.r")

create_marked_synth = function(test_set,
                               remove_id = NULL,
                               remove_id_dat = NULL,
                               k = 10,
                               clock = 'reset') {
  # create k pairs of bootstrap data sets of the test_set and uses them as training
  # sets to create shadow synthesizers of synthetic data. Half of the training
  # data contains remove id and the other half does not. The output is a list, where the
  # first element is a list of the synthetic data sets (of size 2*k),
  # the second element is a
  # vector with its labels (indicating remove_ids presence in the training set)
  # and the final element is remove_id (which can also be randomly sampled)
  
  #k folds
  set_ids = unlist(as.set(test_set$patient_id))
  
  n = length(set_ids)
  
  if (is.null(remove_id)) {
    # draw a random data point to be removed
    x_rem = sample(set_ids, 1)
  }
  else{
    x_rem = remove_id
  }
  
  
  data_rem = test_set[test_set$patient_id != x_rem, ]
  set_ids_rem = setdiff(set_ids, x_rem)
  
  # empty lists for the synthetic data sets and their labels
  boot_syn_list = list()
  label_boot_syn = list()
  
  for (i in 1:k) {
    boot_dat = sample_replacement_pat(data_rem, set_ids_rem, n - 1)
    # add x_rem back in
    boot_dat = rbind(boot_dat, test_set[test_set$patient_id == x_rem, ])
    
    
    # create synthetic data based on boot_dat
    boot_pat = dat_patients(boot_dat)
    boot_syn_pat = gen_patients(boot_pat)
    boot_weibull = fit_weibull(boot_dat, clock)
    boot_seq = gen_sequence_man(boot_syn_pat, boot_weibull, clock = clock)
    boot_syn_dat = seq_to_dataset(boot_seq, boot_dat, boot_syn_pat)
    # 1 to label that x_rem is included
    boot_syn_list = append(boot_syn_list, list(boot_syn_dat))
    label_boot_syn = append(label_boot_syn, 1)
  }
  for (i in 1:k) {
    boot_dat = sample_replacement_pat(data_rem, set_ids_rem, n)
    
    # create synthetic data based on boot_dat_rem
    boot_pat = dat_patients(boot_dat)
    boot_syn_pat = gen_patients(boot_pat)
    boot_weibull = fit_weibull(boot_dat, clock)
    boot_seq = gen_sequence_man(boot_syn_pat, boot_weibull, clock = clock)
    boot_syn_dat = seq_to_dataset(boot_seq, boot_dat, boot_syn_pat)
    # 0 to label that x_rem has been removed
    boot_syn_list = append(boot_syn_list, list(boot_syn_dat))
    label_boot_syn = append(label_boot_syn, 0)
    
  }
  #label_boot_syn = as.list(as.factor(unlist(label_boot_syn)))
  
  return (list(boot_syn_list, label_boot_syn, x_rem))
}


sample_replacement_pat = function(dat, set_ids, n) {
  # sample with replacement from a data set
  boot_ids = sample(set_ids, n,  replace = TRUE)
  boot_dat = dat[dat$patient_id %in% boot_ids,]
  # if patients were drawn more than once, they need to be added again
  freq = as.data.frame(table(boot_ids))
  freq_rep = freq[freq$Freq >= 2,]
  if (nrow(freq_rep) > 0) {
    # number of rows needed to be added
    
    # empty df of approx. dim to be appended to boot_dat
    app_rows = data.frame(matrix(nrow = nrow(dat) * 2 + 20, ncol = ncol(dat)))
    colnames(app_rows) = colnames(dat)
    
    # new ids begin at the newest 1000th
    new_ids_from = ceiling(max(dat$patient_id) / 1000) * 1000
    s = 1
    for (i in 1:nrow(freq_rep)) {
      ith_id = freq_rep[i, 1]
      t = nrow(dat[dat$patient_id == ith_id, ])
      for (j in 1:(freq_rep[i, 2] - 1)) {
        app_rows[s:(s + t - 1), ] = dat[dat$patient_id == ith_id, ]
        
        # new id to distinguish the duplicates
        app_rows[s:(s + t - 1), 1] = app_rows[s:(s + t - 1), 1] + j * new_ids_from
        s = s + t
        if (s > nrow(app_rows)) {
          stop("The data frame is too small.")
        }
        
      }
    }
    boot_dat = rbind.data.frame(boot_dat, app_rows[1:s - 1, ])
    
  }
  
  
  return (boot_dat)
}

combined_df_lab = function(df_list, label_list, split = length(df_list)) {
  # shuffles the data frames and labels returned by create_marked_synth.
  # can also split into test and training
  # adds the label as a column in the data frames
  
  # shuffle the lists
  shuffle = sample(1:length(df_list))
  df_list = df_list[shuffle]
  label_list = label_list[shuffle]
  # find dimensions of the combined data frame
  row_nr = 0
  
  for (i in 1:split) {
    row_nr = row_nr + nrow(df_list[[i]])
  }
  col_nr = ncol(df_list[[1]]) + 1
  
  comb_df_train = as.data.frame(matrix(nrow = row_nr, ncol = col_nr))
  
  row_nr = 1
  for (i in 1:split) {
    row_nr_new = row_nr + nrow(df_list[[i]])
    comb_df_train[row_nr:(row_nr_new - 1), -ncol(comb_df_train)] = df_list[[i]]
    comb_df_train[row_nr:(row_nr_new - 1), ncol(comb_df_train)] = rep(label_list[[i]],
                                                                      nrow(df_list[[i]]))
    row_nr = row_nr_new
  }
  
  colnames(comb_df_train) = c(colnames(df_list[[1]]), 'label')
  if (split < length(df_list)) {
    row_nr = 0
    for (i in (split + 1):length(df_list)) {
      row_nr = row_nr + nrow(df_list[[i]])
    }
    col_nr = ncol(df_list[[1]]) + 1
    
    comb_df_test = as.data.frame(matrix(nrow = row_nr, ncol = col_nr))
    
    row_nr = 1
    for (i in (split + 1):length(df_list)) {
      row_nr_new = row_nr + nrow(df_list[[i]])
      comb_df_test[row_nr:(row_nr_new - 1), -ncol(comb_df_test)] = df_list[[i]]
      comb_df_test[row_nr:(row_nr_new - 1),
                   ncol(comb_df_test)] = rep(label_list[[i]],
                                             nrow(df_list[[i]]))
      row_nr = row_nr_new
    }
    colnames(comb_df_test) = c(colnames(df_list[[1]]), 'label')
    return(list(comb_df_train, comb_df_test))
  }
  else{
    return(list(comb_df_train))
  }
}

train_classifier_mi = function(comb_df_train, comb_df_test = NULL) {
  # trains a random forest classifier, which predicts the label based on
  # the data set
  # NOTE does not consider the data sets as a whole, but row by row
  red_x_train = comb_df_train[, c(1, 4:7, 9:11)]
  if (is.null(comb_df_test)) {
    m = randomForest(y = as.factor(comb_df_train[, 13]),
                     x = red_x_train,
                     data = comb_df_train)
  }
  else{
    red_x_test = comb_df_test[, c(1, 4:7, 9:11)]
    m = randomForest(
      y = as.factor(comb_df_train[, 13]),
      x = red_x_train,
      data = comb_df_train,
      xtest = red_x_test,
      ytest = as.factor(comb_df_test[, 13])
    )
  }
  return (m)
}

predict_synth = function(classifier, synth_dat) {
  pred = mean(as.numeric(predict(classifier, synth_dat))) - 1
  return (pred)
}
