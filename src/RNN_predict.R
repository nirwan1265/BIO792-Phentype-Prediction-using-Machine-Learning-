RNN_func <- function(train_pheno, train_marker, test_marker) {
  # Define the RNN model
  model <- keras_model_sequential()
  model %>% 
    layer_lstm(units = 64, input_shape = c(1, ncol(train_marker))) %>%
    layer_dense(units = 32, activation = "relu") %>%
    layer_dense(units = 1, activation = "linear")
  
  # Compile the model
  model %>% compile(
    loss = "mse",
    optimizer = optimizer_adam(learning_rate = 0.001)
  )
  
  # Reshape the input data
  train_marker <- array_reshape(train_marker, c(nrow(train_marker), 1, ncol(train_marker)))
  test_marker <- array_reshape(test_marker, c(nrow(test_marker), 1, ncol(test_marker)))
  
  # Fit the model to the training data
  history <- model %>% fit(
    x = train_marker,
    y = train_pheno,
    epochs = 50,
    batch_size = 32,
    validation_split = 0.2,
    verbose = 2
  )
  
  # Make predictions on the test and train data
  test_pred <- model %>% predict(test_marker)
  train_pred <- model %>% predict(train_marker)
  return_value = list("test_predicted"=test_pred, "train_predicted"=train_pred, "model"=model)
  # Return the test predictions
  return(return_value)
}

#
#This RNN uses a single LSTM layer with 32 units, followed by a fully connected layer with 64 units 
#and a dropout rate of 0.1. The output layer is a single unit with a linear activation function. 
#The input data is reshaped to have a time dimension of 1, and the batch size is determined based 
#on the size of the training data. 
#The rest of the code is similar to the RNN and DNN functions, with the model being compiled and trained, 
#and predictions being made on the training and test data

# sol_VL
RNN_sol_VL <- RNN_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- RNN_sol_VL$test_pred
train_predicted <- RNN_sol_VL$train_pred


cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

summary(sol_VL_test_phenotype)
summary(test_predicted)

sol_VL_train_test_rnn <- cbind(sol_VL_test_phenotype, test_predicted)




