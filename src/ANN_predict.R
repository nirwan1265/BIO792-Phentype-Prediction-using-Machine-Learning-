ANN_func <- function(train_phenotype, train_marker, test_marker) {
  # Define input shape
  input_shape <- dim(train_marker)[2]
  
  # Define model architecture
  model <- keras_model_sequential() %>%
    layer_dense(units = 128, activation = "relu", input_shape = input_shape) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 1, activation = "sigmoid")
  
  # Compile the model
  model %>% compile(
    loss = "binary_crossentropy",
    optimizer = "adam",
    metrics = c("accuracy")
  )
  
  # Fit the model to the training data
  history <- model %>% fit(
    x = train_marker, y = train_phenotype,
    epochs = 50,
    batch_size = 50,
    shuffle = TRUE,
    validation_split = 0.2,
    verbose = 0
  )
  
  # Use the trained model to predict on the test data
  predicted_phenotype <- predict(model, test_marker)
  
  # Return a list of values
  return(list(predicted_phenotype = predicted_phenotype, history = history))
  # # Use the trained model to predict on the test data
  # train_pred <- model %>% predict(train_marker)
  # test_pred <- model %>% predict(test_marker)
  # return_value = list("test_predicted"=test_pred, "train_predicted"=train_pred, "model"=model
  #                     , history = history)
  # return(return_value)
  
}

# sol_VL
ANN_sol_VL <- ANN_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- ANN_sol_VL$predicted_phenotype
train_predicted <- ANN_sol_VL$train_pred

summary(sol_VL_test_phenotype)
summary(test_predicted)

cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

