RF_func = function(train_pheno, train_geno, test_geno){
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, 
                                        mtry=round(dim(train_geno)[1]/3), ntree=200, 
                                        metric="RMSE", maxnodes=NULL)
  train_predicted <- predict(train_RF, newdata = train_geno)
  val_predicted <- predict(train_RF, newdata = test_geno)
  return_value = list("val_predicted"=val_predicted, "train_predicted"=train_predicted, "model"=train_RF)
  return(return_value)
}
typeof(sol_Mo_test_marker)
class(sol_Mo_test_marker)
dim(sol_Mo_test_marker)
val_predicted <- predict(RF_sol_Mo, newdata = sol_Mo_test_marker)

# RF_func = function (train_pheno, train_geno, test_geno, check){
#   train_geno = (train_geno + 1) / 3
#   test_geno = (test_geno + 1) / 3
#   if (check == "FM"){
#     train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, mtry=round(dim(train_geno)[1]/3), ntree=200, metric="RMSE", maxnodes=NULL)
#   } else {
#     train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, mtry=round(dim(train_geno)[1]), ntree=100, metric="RMSE", maxnodes=NULL)
#   }
#   train_predicted <- predict(train_RF, newdata = train_geno)
#   val_predicted <- predict(train_RF, newdata = test_geno)
#   return_value = list("val_predicted"=val_predicted, "train_predicted"=train_predicted, "model"=train_RF)
#   return(return_value)
# }

# sol_Mo
RF_sol_Mo <- RF_func(sol_Mo_train_phenotype, sol_Mo_train_marker, sol_Mo_test_marker)

test_predicted <- RF_sol_Mo$val_predicted
train_predicted <- RF_sol_Mo$train_predicted

cor(sol_Mo_train_phenotype,train_predicted)
cor(sol_Mo_test_phenotype, test_predicted)

summary(test_predicted)
summary(train_predicted)



sol_Mo_train_test_RF <- cbind(sol_Mo_test_phenotype, test_predicted)


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


