SVM_func <- function(train_pheno, train_geno, test_geno) {
  
  # Normalize genotype data to [0, 1] range
  train_geno_norm <- apply(train_geno, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  test_geno_norm <- apply(test_geno, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  
  # Train SVM model with radial basis function kernel
  svm_model <- svm(train_pheno ~., data = data.frame(train_geno_norm), kernel = "radial")
  
  # Make predictions on training and testing data
  train_predicted <- predict(svm_model, data.frame(train_geno_norm))
  val_predicted <- predict(svm_model, data.frame(test_geno_norm))
  
  return_value <- list("val_predicted" = val_predicted, "train_predicted" = train_predicted, "model" = svm_model)
  return(return_value)
}

# stp
SVM_stp <- SVM_func(stp_train_phenotype, stp_train_marker, stp_test_marker)

test_predicted <- SVM_stp$val_predicted
train_predicted <- SVM_stp$train_predicted

summary(test_predicted)
summary(stp_test_phenotype)


cor(stp_train_phenotype,train_predicted)
cor(stp_test_phenotype, test_predicted)

stp_train_test_SVM <- cbind(stp_test_phenotype, test_predicted)
