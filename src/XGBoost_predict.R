xgb_func <- function(train_pheno, train_geno, test_geno, nrounds = 1000, eta = 0.01, max_depth = 3) {
  # convert train_geno and test_geno to DMatrix objects
  train_geno <- xgb.DMatrix(train_geno, label = train_pheno)
  test_geno <- xgb.DMatrix(test_geno)
  
  # Define the XGBoost model
  xgb_model <- xgboost(
    data = train_geno,
    nrounds = nrounds,
    objective = "reg:squarederror",
    eta = eta,
    max_depth = max_depth,
    verbose = 0
  )
  
  # Predict the phenotypes for training and test data
  train_predicted <- predict(xgb_model, train_geno)
  val_predicted <- predict(xgb_model, test_geno)
  
  # Return predicted values and model object
  return_value <- list("train_predicted" = train_predicted, "test_predicted" = test_predicted, "model" = xgb_model)
  return(return_value)
}
test_predicted <- predict(xgb_sol_Mo, as.matrix(sol_Mo_test_marker))

class(sol_Mo_test_marker)
# sol_VL
xgb_sol_VL <- xgb_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- xgb_sol_VL$test_predicted
train_predicted <- xgb_sol_VL$train_predicted

summary(test_predicted)
summary(sol_VL_test_phenotype)


cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)
