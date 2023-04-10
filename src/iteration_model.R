# Set the number of iterations to run the loop
num_iter <- 100

# Correlation
train_corr_sol_VL <- vector("numeric", num_iter)
test_corr_sol_VL <- vector("numeric", num_iter)

train_corr_tot <- vector("numeric", num_iter)
test_corr_tot <- vector("numeric", num_iter)

train_corr_lab <- vector("numeric", num_iter)
test_corr_lab <- vector("numeric", num_iter)


# Create empty vectors to store the predicted values for each iteration
train_preds_sol_VL <- matrix(NA, nrow = nrow(sol_VL), ncol = num_iter)
test_preds_sol_VL <- matrix(NA, nrow = nrow(sol_VL), ncol = num_iter)

train_preds_tot <- matrix(NA, nrow = nrow(tot), ncol = num_iter)
test_preds_tot <- matrix(NA, nrow = nrow(tot), ncol = num_iter)

train_preds_lab <- matrix(NA, nrow = nrow(lab), ncol = num_iter)
test_preds_lab <- matrix(NA, nrow = nrow(lab), ncol = num_iter)


# Run the loop
for (i in 1:num_iter) {
  # Randomly sample markers and phenotypes for the training set
  tot_train <- as.matrix(sample(1:nrow(tot_markers), 1035))
  tot_test <- setdiff(1:nrow(tot_markers), tot_train)
  tot_train_phenotype <- as.matrix(tot[tot_train, ])
  tot_train_marker <- as.matrix(tot_markers[tot_train, ], K = NULL)
  # Removing NAs from training and putting them into testing
  nas <- which(is.na(tot_train_phenotype))
  tot_train_phenotype <- as.matrix(tot_train_phenotype[-nas, ])
  tot_train_marker <- as.matrix(tot_train_marker[-nas, ])
  tot_train <- tot_train[-nas,]
  tot_test <- sort(c(tot_test, nas))
  tot_test_phenotype <- as.matrix(tot[tot_test, ])
  tot_test_marker <- as.matrix(tot_markers[tot_test, ], K = NULL)
  
  
  lab_train <- as.matrix(sample(1:nrow(lab_markers), 1035))
  lab_test <- setdiff(1:nrow(lab_markers), lab_train)
  lab_train_phenotype <- as.matrix(lab[lab_train, ])
  lab_train_marker <- as.matrix(lab_markers[lab_train, ], K = NULL)
  # Removing NAs from training and putting them into testing
  nas <- which(is.na(lab_train_phenotype))
  lab_train_phenotype <- as.matrix(lab_train_phenotype[-nas, ])
  lab_train_marker <- as.matrix(lab_train_marker[-nas, ])
  lab_train <- lab_train[-nas,]
  lab_test <- sort(c(lab_test, nas))
  lab_test_phenotype <- as.matrix(lab[lab_test, ])
  lab_test_marker <- as.matrix(lab_markers[lab_test, ], K = NULL)

  sol_VL_train_phenotype <- as.matrix(sol_VL[sol_VL_train, ])
  sol_VL_train_marker <- as.matrix(sol_VL_markers[sol_VL_train, ], K = NULL)
  sol_VL_test_phenotype <- as.matrix(sol_VL[sol_VL_test, ])
  sol_VL_test_marker <- as.matrix(sol_VL_markers[sol_VL_test, ], K = NULL)
  
    # Fit the rrBLUP model to the training data
  # RRb_tot <- RRb_func(tot_train_phenotype, tot_train_marker, tot_test_marker)
  # RRb_sol_VL <- RRb_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)
  # RRb_lab <- RRb_func(lab_train_phenotype, lab_train_marker, lab_test_marker)
  #RRb_tot <- DNN_func(tot_train_phenotype, tot_train_marker, tot_test_marker)
  RRb_sol_VL <- RNN_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)
  #RRb_lab <- DNN_func(lab_train_phenotype, lab_train_marker, lab_test_marker)
  

  # Extract the predicted values for the training and test sets
  train_preds_sol_VL[1:nrow(RRb_sol_VL$train_predicted), i] <- RRb_sol_VL$train_predicted
  test_preds_sol_VL[1:nrow(RRb_sol_VL$val_predicted), i] <- RRb_sol_VL$val_predicted
  # Calculate correlation without NA values
  train_corr_sol_VL[i] <- cor(sol_VL_train_phenotype, train_preds_sol_VL[1:nrow(RRb_sol_VL$train_predicted), i], use = "pairwise.complete.obs")
  test_corr_sol_VL[i] <- cor(sol_VL_test_phenotype, test_preds_sol_VL[1:nrow(RRb_sol_VL$val_predicted), i], use = "pairwise.complete.obs")
   
  # train_preds_tot[1:nrow(RRb_tot$train_predicted), i] <- RRb_tot$train_predicted
  # test_preds_tot[1:nrow(RRb_tot$val_predicted), i] <- RRb_tot$val_predicted
  # Calculate correlation without NA values
  # train_corr_tot[i] <- cor(tot_train_phenotype, train_preds_tot[1:nrow(RRb_tot$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_tot[i] <- cor(tot_test_phenotype, test_preds_tot[1:nrow(RRb_tot$val_predicted), i], use = "pairwise.complete.obs")
  # 
  # train_preds_lab[1:nrow(RRb_lab$train_predicted), i] <- RRb_lab$train_predicted
  # test_preds_lab[1:nrow(RRb_lab$val_predicted), i] <- RRb_lab$val_predicted
  # Calculate correlation without NA values
  # train_corr_lab[i] <- cor(lab_train_phenotype, train_preds_lab[1:nrow(RRb_lab$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_lab[i] <- cor(lab_test_phenotype, test_preds_lab[1:nrow(RRb_lab$val_predicted), i], use = "pairwise.complete.obs")
}

# Calculate mean correlation
mean_train_corr_sol_VL <- mean(train_corr_sol_VL)
mean_test_corr_sol_VL <- mean(test_corr_sol_VL)

mean_train_corr_tot <- mean(train_corr_tot)
mean_test_corr_tot <- mean(test_corr_tot)

mean_train_corr_lab <- mean(train_corr_lab)
mean_test_corr_lab <- mean(test_corr_lab)
