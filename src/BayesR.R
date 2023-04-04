BayesR_func_lasso <- function(train_pheno, train_geno, test_geno, alpha = 1) {
  
  # Center and scale the genotype matrix
  train_geno <- scale(train_geno, center = TRUE, scale = TRUE)
  test_geno <- scale(test_geno, center = TRUE, scale = TRUE)
  
  # Fit LASSO model to select optimal lambda
  cvfit <- cv.glmnet(train_geno, train_pheno, alpha = alpha, family = "gaussian", standardize = FALSE, intercept = FALSE)
  lambda <- cvfit$lambda.min
  
  # Define prior for model weights
  prior_w <- rep(0, ncol(train_geno) + 1)
  prior_w[1] <- mean(train_pheno)
  
  # Compute the BayesR model
  library(BayesR)
  br <- BayesR(y = train_pheno, X = train_geno, S = 1, prior_w = prior_w, lambda = lambda, alpha = alpha)
  
  # Predict phenotypes for training and test data
  train_pred <- predict(br, train_geno)
  test_pred <- predict(br, test_geno)
  
  # Return predicted values and model object
  return_value <- list("train_predicted" = train_pred, "test_predicted" = test_pred, "model" = br)
  return(return_value)
}

# sol_VL
BayesR_sol_VL <- BayesR_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- BayesR_sol_VL$test_predicted
train_predicted <- BayesR_sol_VL$train_predicted

summary(test_predicted)
summary(sol_VL_test_phenotype)


cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

sol_VL_train_test_BayesR <- cbind(sol_VL_test_phenotype, test_predicted)
