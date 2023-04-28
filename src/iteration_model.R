# Set the number of iterations to run the loop
num_iter <- 25

# Correlation
train_corr_stp <- vector("numeric", num_iter)
test_corr_stp <- vector("numeric", num_iter)

train_corr_stp <- vector("numeric", num_iter)
test_corr_stp <- vector("numeric", num_iter)

train_corr_lab <- vector("numeric", num_iter)
test_corr_lab <- vector("numeric", num_iter)


# Create empty vectors to store the predicted values for each iteration
train_preds_stp <- matrix(NA, nrow = nrow(stp), ncol = num_iter)
test_preds_stp <- matrix(NA, nrow = nrow(stp), ncol = num_iter)

train_preds_stp <- matrix(NA, nrow = nrow(stp), ncol = num_iter)
test_preds_stp <- matrix(NA, nrow = nrow(stp), ncol = num_iter)

train_preds_lab <- matrix(NA, nrow = nrow(lab), ncol = num_iter)
test_preds_lab <- matrix(NA, nrow = nrow(lab), ncol = num_iter)


# Initialize lists to store test predictions
test_preds_stp_list <- vector("list", num_iter)
test_preds_PBR_list <- vector("list", num_iter)
test_preds_PNZ_list <- vector("list", num_iter)
test_preds_sol_Mo_list <- vector("list", num_iter)




# Run the loop
for (i in 1:num_iter) {
  # Randomly sample markers and phenotypes for the training set
  stp_train <- as.matrix(sample(1:nrow(stp_markers), 1035))
  stp_test <- setdiff(1:nrow(stp_markers), stp_train)
  stp_train_phenotype <- as.matrix(stp[stp_train, ])
  stp_train_marker <- as.matrix(stp_markers[stp_train, ], K = NULL)
  # Removing NAs from training and putting them into testing
  nas <- which(is.na(stp_train_phenotype))
  stp_train_phenotype <- as.matrix(stp_train_phenotype[-nas, ])
  stp_train_marker <- as.matrix(stp_train_marker[-nas, ])
  stp_train <- stp_train[-nas,]
  stp_test <- sort(c(stp_test, nas))
  stp_test_phenotype <- as.matrix(stp[stp_test, ])
  stp_test_marker <- as.matrix(stp_markers[stp_test, ], K = NULL)
  
  
  # lab_train <- as.matrix(sample(1:nrow(lab_markers), 1035))
  # lab_test <- setdiff(1:nrow(lab_markers), lab_train)
  # lab_train_phenotype <- as.matrix(lab[lab_train, ])
  # lab_train_marker <- as.matrix(lab_markers[lab_train, ], K = NULL)
  # Removing NAs from training and putting them into testing
  # nas <- which(is.na(lab_train_phenotype))
  # lab_train_phenotype <- as.matrix(lab_train_phenotype[-nas, ])
  # lab_train_marker <- as.matrix(lab_train_marker[-nas, ])
  # lab_train <- lab_train[-nas,]
  # lab_test <- sort(c(lab_test, nas))
  # lab_test_phenotype <- as.matrix(lab[lab_test, ])
  # lab_test_marker <- as.matrix(lab_markers[lab_test, ], K = NULL)
  # 
  # stp_train_phenotype <- as.matrix(stp[stp_train, ])
  # stp_train_marker <- as.matrix(stp_markers[stp_train, ], K = NULL)
  # stp_test_phenotype <- as.matrix(stp[stp_test, ])
  # stp_test_marker <- as.matrix(stp_markers[stp_test, ], K = NULL)
  # 
  #   # Fit the RRbLUP model to the training data
  # RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  # RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  # RRb_lab <- RRb_func(lab_train_phenotype, lab_train_marker, lab_test_marker)
    RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  # RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  #RRb_lab <- RRb_func(lab_train_phenotype, lab_train_marker, lab_test_marker)
  

  # Extract the predicted values for the training and test sets
  # train_preds_stp[1:nrow(RRb_stp$train_predicted), i] <- RRb_stp$train_predicted
  # test_preds_stp[1:nrow(RRb_stp$val_predicted), i] <- RRb_stp$val_predicted
  # # Calculate correlation without NA values
  # train_corr_stp[i] <- cor(stp_train_phenotype, train_preds_stp[1:nrow(RRb_stp$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_stp[i] <- cor(stp_test_phenotype, test_preds_stp[1:nrow(RRb_stp$val_predicted), i], use = "pairwise.complete.obs")
  #  
  train_preds_stp[1:nrow(RRb_stp$train_predicted), i] <- RRb_stp$train_predicted
  test_preds_stp[1:nrow(RRb_stp$val_predicted), i] <- RRb_stp$val_predicted
  #Calculate correlation without NA values
  train_corr_stp[i] <- cor(stp_train_phenotype, train_preds_stp[1:nrow(RRb_stp$train_predicted), i], use = "pairwise.complete.obs")
  test_corr_stp[i] <- cor(stp_test_phenotype, test_preds_stp[1:nrow(RRb_stp$val_predicted), i], use = "pairwise.complete.obs")

  # train_preds_lab[1:nrow(RRb_lab$train_predicted), i] <- RRb_lab$train_predicted
  # test_preds_lab[1:nrow(RRb_lab$val_predicted), i] <- RRb_lab$val_predicted
  # Calculate correlation without NA values
  # train_corr_lab[i] <- cor(lab_train_phenotype, train_preds_lab[1:nrow(RRb_lab$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_lab[i] <- cor(lab_test_phenotype, test_preds_lab[1:nrow(RRb_lab$val_predicted), i], use = "pairwise.complete.obs")
  
  # Initialize empty matrices to store the test predictions for each phenotype
  stp_pred_matrix <- matrix(NA, nrow = nrow(stp), 1)
  # PBR_pred_matrix <- matrix(NA, nrow = nrow(PBR), 1)
  # PNZ_pred_matrix <- matrix(NA, nrow = nrow(PNZ), 1)
  # sol_Mo_pred_matrix <- matrix(NA, nrow = nrow(sol_Mo), 1)
  # 
  # Fill the matrices with the test predictions based on the test indices
  stp_pred_matrix[stp_test, 1] <- RRb_stp$val_pred
  # PBR_pred_matrix[PBR_test, 1] <- RRb_PBR$val_predicted
  # PNZ_pred_matrix[PNZ_test, 1] <- RRb_PNZ$val_predicted
  # sol_Mo_pred_matrix[sol_Mo_test, 1] <- RRb_sol_Mo$val_predicted
  # 
  # Store the filled matrices in the respective lists
  test_preds_stp_list[[i]] <- stp_pred_matrix
  
}

mean_train_corr_stp <- mean(train_corr_stp)
mean_test_corr_stp <- mean(test_corr_stp)

mean_train_corr_stp <- mean(train_corr_stp)
mean_test_corr_stp <- mean(test_corr_stp)



for (i in 1:num_iter) {
  # Randomly sample markers and phenotypes for the training set
  stp_train <- as.matrix(sample(1:nrow(stp_markers), 1035))
  stp_test <- setdiff(1:nrow(stp_markers), stp_train)
  stp_train_phenotype <- as.matrix(stp[stp_train, ])
  stp_train_marker <- as.matrix(stp_markers[stp_train, ], K = NULL)
  # Removing NAs from training and putting them into testing
  nas <- which(is.na(stp_train_phenotype))
  stp_train_phenotype <- as.matrix(stp_train_phenotype[-nas, ])
  stp_train_marker <- as.matrix(stp_train_marker[-nas, ])
  stp_train <- stp_train[-nas,]
  stp_test <- sort(c(stp_test, nas))
  stp_test_phenotype <- as.matrix(stp[stp_test, ])
  stp_test_marker <- as.matrix(stp_markers[stp_test, ], K = NULL)
  
  
  # lab_train <- as.matrix(sample(1:nrow(lab_markers), 1035))
  # lab_test <- setdiff(1:nrow(lab_markers), lab_train)
  # lab_train_phenotype <- as.matrix(lab[lab_train, ])
  # lab_train_marker <- as.matrix(lab_markers[lab_train, ], K = NULL)
  # Removing NAs from training and putting them into testing
  # nas <- which(is.na(lab_train_phenotype))
  # lab_train_phenotype <- as.matrix(lab_train_phenotype[-nas, ])
  # lab_train_marker <- as.matrix(lab_train_marker[-nas, ])
  # lab_train <- lab_train[-nas,]
  # lab_test <- sort(c(lab_test, nas))
  # lab_test_phenotype <- as.matrix(lab[lab_test, ])
  # lab_test_marker <- as.matrix(lab_markers[lab_test, ], K = NULL)
  # 
  # stp_train_phenotype <- as.matrix(stp[stp_train, ])
  # stp_train_marker <- as.matrix(stp_markers[stp_train, ], K = NULL)
  # stp_test_phenotype <- as.matrix(stp[stp_test, ])
  # stp_test_marker <- as.matrix(stp_markers[stp_test, ], K = NULL)
  # 
  #   # Fit the RRbLUP model to the training data
  # RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  # RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  # RRb_lab <- RRb_func(lab_train_phenotype, lab_train_marker, lab_test_marker)
  #RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  # RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  #RRb_lab <- RRb_func(lab_train_phenotype, lab_train_marker, lab_test_marker)
  
  # Wrap the RRb_func() call inside a tryCatch() block
  # Wrap the RRb_func() call inside a tryCatch() block
  error_flag <- FALSE
  tryCatch({
    RRb_stp <- RRb_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  }, error = function(e) {
    error_flag <- TRUE
  })
  
  # If an error occurred during RRb_func(), skip the rest of the current iteration
  if (error_flag) {
    next
  }
  
  
  
  # Extract the predicted values for the training and test sets
  # train_preds_stp[1:nrow(RRb_stp$train_predicted), i] <- RRb_stp$train_predicted
  # test_preds_stp[1:nrow(RRb_stp$val_predicted), i] <- RRb_stp$val_predicted
  # # Calculate correlation without NA values
  # train_corr_stp[i] <- cor(stp_train_phenotype, train_preds_stp[1:nrow(RRb_stp$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_stp[i] <- cor(stp_test_phenotype, test_preds_stp[1:nrow(RRb_stp$val_predicted), i], use = "pairwise.complete.obs")
  #  
  train_preds_stp[1:nrow(RRb_stp$train_predicted), i] <- RRb_stp$train_predicted
  test_preds_stp[1:nrow(RRb_stp$val_predicted), i] <- RRb_stp$val_predicted
  #Calculate correlation without NA values
  #train_corr_stp[i] <- cor(stp_train_phenotype, train_preds_stp[1:nrow(RRb_stp$train_predicted), i], use = "pairwise.complete.obs")
  #test_corr_stp[i] <- cor(stp_test_phenotype, test_preds_stp[1:nrow(RRb_stp$val_predicted), i], use = "pairwise.complete.obs")
  train_corr_stp[i] <- cor(stp_train_phenotype[1:nrow(RRb_stp$train_predicted), ], train_preds_stp[1:nrow(RRb_stp$train_predicted), i], use = "pairwise.complete.obs")
  test_corr_stp[i] <- cor(stp_test_phenotype[1:nrow(RRb_stp$val_predicted), ], test_preds_stp[1:nrow(RRb_stp$val_predicted), i], use = "pairwise.complete.obs")
  
    # train_preds_lab[1:nrow(RRb_lab$train_predicted), i] <- RRb_lab$train_predicted
  # test_preds_lab[1:nrow(RRb_lab$val_predicted), i] <- RRb_lab$val_predicted
  # Calculate correlation without NA values
  # train_corr_lab[i] <- cor(lab_train_phenotype, train_preds_lab[1:nrow(RRb_lab$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_lab[i] <- cor(lab_test_phenotype, test_preds_lab[1:nrow(RRb_lab$val_predicted), i], use = "pairwise.complete.obs")
  
  # Initialize empty matrices to store the test predictions for each phenotype
  stp_pred_matrix <- matrix(NA, nrow = nrow(stp), 1)
  # PBR_pred_matrix <- matrix(NA, nrow = nrow(PBR), 1)
  # PNZ_pred_matrix <- matrix(NA, nrow = nrow(PNZ), 1)
  # sol_Mo_pred_matrix <- matrix(NA, nrow = nrow(sol_Mo), 1)
  # 
  # Fill the matrices with the test predictions based on the test indices
  stp_pred_matrix[stp_test, 1] <- RRb_stp$val_pred
  # PBR_pred_matrix[PBR_test, 1] <- RRb_PBR$val_predicted
  # PNZ_pred_matrix[PNZ_test, 1] <- RRb_PNZ$val_predicted
  # sol_Mo_pred_matrix[sol_Mo_test, 1] <- RRb_sol_Mo$val_predicted
  # 
  # Store the filled matrices in the respective lists
  test_preds_stp_list[[i]] <- stp_pred_matrix
  
}
