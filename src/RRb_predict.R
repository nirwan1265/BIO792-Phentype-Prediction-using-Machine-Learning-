RRb_func = function (train_pheno, train_geno, test_geno) {
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  train_BLUP <- mixed.solve(y = train_pheno, Z = train_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)
  train_e = as.matrix(train_BLUP$u)
  train_valid = train_geno %*% train_e
  train_pred = train_valid + c(train_BLUP$beta)
  val_valid = test_geno %*% train_e
  val_pred <- val_valid + c(train_BLUP$beta)
  return_value = list("val_predicted"=val_pred, "train_predicted"=train_pred, "model"=train_BLUP, "u"=train_BLUP$u, "beta"=train_BLUP$beta)
  return(return_value)
}

# sol_VL
RRb_sol_VL <- RRb_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- RRb_sol_VL$val_predicted
train_predicted <- RRb_sol_VL$train_predicted

cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

sol_VL_train_test_RRb <- cbind(sol_VL_test_phenotype, test_predicted)
