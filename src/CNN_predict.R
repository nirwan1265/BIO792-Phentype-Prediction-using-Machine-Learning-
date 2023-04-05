CNN_func = function (train_pheno, train_geno, test_geno) {
  mnum = dim(train_geno)[2]
  train_geno_mnum <- train_geno + 1
  dim(train_geno_mnum) <- c(dim(train_geno_mnum)[1], 1, dim(train_geno_mnum)[2])
  train_geno_mnum <- tensorflow::tf$cast(train_geno_mnum, tf$float32) /3
  train_geno_mnum <- tensorflow::tf$expand_dims(train_geno_mnum, axis=-1L)
  test_geno_mnum = test_geno + 1
  dim(test_geno_mnum) <- c(dim(test_geno_mnum)[1], 1, dim(test_geno_mnum)[2])
  test_geno_mnum <- tensorflow::tf$cast(test_geno_mnum, tf$float32) /3
  test_geno_mnum <- tensorflow::tf$expand_dims(test_geno_mnum, axis=-1L)
  batchs = round(dim(train_geno_mnum)[1]/20)
  if ((batchs %% 2) == 1){batchs = batchs + 1}
  model <- keras_model_sequential() %>%
    layer_conv_2d(filters = 32, kernel_size = c(1,14), strides = c(1,4), padding="same", activation= "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.1, l2 = 0.1), input_shape = shape(1, mnum, 1)) %>%
    layer_dropout(rate = 0.2) %>%
    layer_conv_2d(filters = 16, kernel_size = c(1,10), strides = c(1,3), padding="same", activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_dropout(rate = 0.1) %>%
    layer_conv_2d(filters = 8, kernel_size = c(1,8), strides = c(1,2), padding="same", activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_max_pooling_2d(pool_size = c(1,2))  %>%
    layer_batch_normalization(batch_size = batchs)
  model <- model %>%
    layer_flatten() %>%
    layer_dense(units = 64, activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_dense(units = 32, activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_dense(units = 16, activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_batch_normalization(batch_size = batchs)%>%
    layer_dense(units = 1, activation = "linear")
  model <- model %>% compile(loss = "mse", optimizer = tf$keras$optimizers$legacy$Adamax(learning_rate=0.003, decay = 0.0003), metrics = c("mean_absolute_error", "mean_squared_error"))
  tensorflow::tf$random$set_seed(200)
  tf$compat$v1$set_random_seed(200)
  set.seed(200)
  tensorflow::tf$random$set_seed(200)
  history <- model %>%
    fit(train_geno_mnum, train_pheno, epochs = 500, batch_size = batchs, validation_split = 0.2, verbose = 0,
        callbacks = list(callback_early_stopping(patience = 30), callback_reduce_lr_on_plateau(factor = 0.1)))
  train_pred <- model %>% predict(train_geno_mnum, batch_size = batchs)
  val_pred <- model %>% predict(test_geno_mnum, batch_size = batchs)
  return_value = list("val_predicted"=val_pred, "train_predicted"=train_pred, "model"=model,  "batch"=batchs)
  return(return_value)
}

# sol_VL
CNN_sol_VL <- CNN_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- CNN_sol_VL$val_predicted
train_predicted <- CNN_sol_VL$train_predicted

cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

summary(sol_VL_test_phenotype)
summary(test_predicted)



sol_VL_train_test_CNN <- cbind(sol_VL_test_phenotype, test_predicted)

#dim(SNP_training_data_nz)#[1:5,1:5]

# lab

CNN_lab <- CNN_func(lab_train_phenotype, lab_train_marker, lab_test_marker)

test_predicted <- CNN_lab$val_predicted
train_predicted <- CNN_lab$train_predicted

cor(lab_train_phenotype,train_predicted)
quartz()
cor(lab_test_phenotype, test_predicted)
plot(lab_test_phenotype, test_predicted)



CNN_tot <- CNN_func(tot_train_phenotype, tot_train_marker, tot_test_marker)

test_predicted <- CNN_tot$val_predicted
train_predicted <- CNN_tot$train_predicted

summary(tot_test_phenotype)
summary(test_predicted)

cor(tot_train_phenotype,train_predicted)
quartz()
cor(tot_test_phenotype, test_predicted)
plot(tot_test_phenotype, test_predicted)
