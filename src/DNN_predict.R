DNN_func = function (train_pheno, train_geno, test_geno) {
  train_geno_mnum = train_geno + 1
  train_geno_mnum = tf$cast(train_geno_mnum, tf$float32) / 3
  test_geno_mnum = test_geno + 1
  test_geno_mnum <- tf$cast(test_geno_mnum, tf$float32) /3
  batchs = round(dim(train_geno_mnum)[1]/20)
  if ((batchs %% 2) == 1){batchs = batchs + 1}
  model <- keras_model_sequential()
  model <- model %>%
    layer_dense(units = 256, activation = "relu", input_shape = dim(train_geno_mnum)[2], kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.3, l2 = 0.3)) %>%
    layer_dropout(rate = 0.03) %>%
    layer_dense(units = 128, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.2, l2 = 0.2)) %>%
    layer_dropout(rate = 0.02) %>%
    layer_dense(units = 64, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_dropout(rate = 0.01) %>%
    layer_dense(units = 32, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_batch_normalization(batch_size = batchs) %>%
    layer_dense(units = 16, activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_dense(units = 1, activation = "linear")
  model <- model %>% compile(loss = "mse", optimizer = tf$keras$optimizers$legacy$Adamax(learning_rate=0.001, decay = 0.0003), metrics = c("mean_absolute_error", "mean_squared_error"))
  #model %>% compile(loss = "mse", optimizer = optimizer_adamax(lr=0.001, decay = 0.0003), metrics = c(metric_r2_score, metric_cor))
  tensorflow::tf$random$set_seed(200)
  tf$compat$v1$set_random_seed(200)
  set.seed(200)
  tensorflow::tf$random$set_seed(200)
  history <- model %>%  
    fit(train_geno_mnum, train_pheno, epochs = 250, batch_size = 10, validation_split = 0.2, verbose = 0,
        callbacks = list(callback_early_stopping(patience = 30), callback_reduce_lr_on_plateau(factor = 0.1)))
  train_pred <- model %>% predict(train_geno_mnum, batch_size = batchs)
  val_pred <- model %>% predict(test_geno_mnum, batch_size = batchs)
  return_value = list("val_predicted"=val_pred, "train_predicted"=train_pred, "model"=model, "batch"=batchs)
  return(return_value)
}
# Changed epochs to 100 and batch_size to 10 from batchs
# Not a lot changed from epoch 100 to 1000


# sol_VL
DNN_sol_VL <- DNN_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- DNN_sol_VL$val_predicted
train_predicted <- DNN_sol_VL$train_predicted

summary(sol_VL_test_phenotype)
summary(test_predicted)

cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

sol_VL_train_test_DNN <- cbind(sol_VL_test_phenotype, test_predicted)




# lab
DNN_lab <- DNN_func(lab_train_phenotype, lab_train_marker, lab_test_marker)
test_predicted <- DNN_lab$val_predicted
train_predicted <- DNN_lab$train_predicted

cor(lab_train_phenotype,train_predicted)



# Get row indices that have complete data in both matrices
complete_rows <- complete.cases(lab_test_phenotype, test_predicted)

# Subset matrices to include only complete rows
lab_test_phenotype_complete <- lab_test_phenotype[complete_rows, ]
test_predicted_complete <- test_predicted[complete_rows, ]

# Compute correlation
cor(lab_test_phenotype_complete, test_predicted_complete)





# DNN
DNN_tot <- DNN_func(tot_train_phenotype, tot_train_marker, tot_test_marker)

test_predicted <- DNN_tot$val_predicted
train_predicted <- DNN_tot$train_predicted

summary(tot_test_phenotype)
summary(test_predicted)
