RF_func = function (train_pheno, train_geno, test_geno, check){
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  set.seed(200)
  if (check == "FM"){
    train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, mtry=round(dim(train_geno)[1]/3), ntree=5000, metric="RMSE", maxnodes=NULL)
  } else {
    train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, mtry=round(dim(train_geno)[1]), ntree=100, metric="RMSE", maxnodes=NULL)
  }
  train_pred <- predict(train_RF, newdata = train_geno)
  val_pred <- predict(train_RF, newdata = test_geno)
  return_value = list("val_predicted"=val_pred, "train_predicted"=train_pred, "model"=train_RF)
  return(return_value)
}

# sol_VL
RF_sol_VL <- RF_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker, "FM")

test_predicted <- RF_sol_VL$val_predicted
train_predicted <- RF_sol_VL$train_predicted

cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

sol_VL_train_test_RF <- cbind(sol_VL_test_phenotype, test_predicted)


RF_lab <- RF_func(lab_train_phenotype, lab_train_marker, lab_test_marker, "FM")
test_predicted <- RF_lab$val_predicted
train_predicted <- RF_lab$train_predicted

cor(lab_train_phenotype,train_predicted)
# Get row indices that have complete data in both matrices
complete_rows <- complete.cases(lab_test_phenotype, test_predicted)

# Subset matrices to include only complete rows
lab_test_phenotype_complete <- lab_test_phenotype[complete_rows, ]
test_predicted_complete <- test_predicted[complete_rows, ]

# Compute correlation
cor(lab_test_phenotype, test_predicted, use = "complete.obs")


RF_tot <- RF_func(tot_train_phenotype, tot_train_marker, tot_test_marker, "FM")
test_predicted <- RF_tot$val_predicted
train_predicted <- RF_tot$train_predicted

summary(tot_test_phenotype)
summary(test_predicted)

cor(tot_train_phenotype,train_predicted)
cor(tot_test_phenotype, test_predicted, use = "complete.obs")


