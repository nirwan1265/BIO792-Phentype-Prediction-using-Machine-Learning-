SVM_func <- function(train_pheno, train_geno, test_geno) {
  
  # Normalize genotype data to [0, 1] range
  train_geno_norm <- apply(train_geno, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  test_geno_norm <- apply(test_geno, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  
  # Train SVM model with radial basis function kernel
  svm_model <- svm(train_pheno ~., data = data.frame(train_geno_norm), kernel = "radial")
  
  # Make predictions on training and testing data
  train_pred <- predict(svm_model, data.frame(train_geno_norm))
  test_pred <- predict(svm_model, data.frame(test_geno_norm))
  
  return_value <- list("val_predicted" = test_pred, "train_predicted" = train_pred, "model" = svm_model)
  return(return_value)
}

# sol_VL
SVM_sol_VL <- SVM_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- SVM_sol_VL$val_predicted
train_predicted <- SVM_sol_VL$train_predicted

summary(test_predicted)
summary(sol_VL_test_phenotype)


cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

sol_VL_train_test_SVM <- cbind(sol_VL_test_phenotype, test_predicted)
