autoencoder_func <- function(train_data, test_data, encoding_dim = 64, epochs = 50) {
  # Define the input shape
  input_shape <- dim(train_data)[2]
  
  # Define the encoder network
  encoder <- keras_model_sequential() %>%
    layer_dense(units = 128, activation = "relu", input_shape = input_shape) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = as.integer(encoding_dim), activation = "relu") 
  
  # Define the decoder network
  decoder <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = "relu", input_shape = encoding_dim) %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = input_shape, activation = "sigmoid")
  
  # Define the autoencoder model
  autoencoder <- keras_model_sequential() %>%
    add(encoder) %>%
    add(decoder)
  
  # Compile the model
  autoencoder %>% compile(
    optimizer = "adam",
    loss = "mse"
  )
  
  # Train the model
  autoencoder %>% fit(
    train_data, train_data,
    epochs = epochs,
    batch_size = 32,
    shuffle = TRUE,
    validation_data = list(test_data, test_data)
  )
  
  # Extract the encoded data
  encoded_data <- encoder %>% predict(train_data)
  
  # Extract the decoded data
  decoded_data <- decoder %>% predict(encoded_data)
  
  # Return the encoded and decoded data as a list
  return(list(encoded_data = encoded_data, decoded_data = decoded_data))
}

autoencoder_func <- function(train_phenotype, train_marker, test_marker) {
  # Define input shape
  input_shape <- dim(train_marker)[2]
  
  # Define encoding dimension
  encoding_dim <- 10
  
  # Define model architecture
  model <- keras_model_sequential() %>%
    layer_dense(units = 128, activation = "relu", input_shape = input_shape) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = encoding_dim, activation = "relu") %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = input_shape, activation = "sigmoid")
  
  # Compile the model
  model %>% compile(
    loss = "mean_squared_error",
    optimizer = "adam"
  )
  
  # Fit the model to the training data
  history <- model %>% fit(
    x = train_marker, y = train_marker,
    epochs = 50,
    batch_size = 128,
    shuffle = TRUE,
    validation_split = 0.2,
    verbose = 0
  )
  
  # Use the trained model to predict on the test data
  predicted_marker <- predict(model, test_marker)
  
  # Compute the mean squared error between the predicted and actual marker data
  mse <- mean((predicted_marker - test_marker)^2)
  
  # Return a list of values
  return(list(predicted_marker = predicted_marker, mse = mse, history = history))
}



autoencoder_func <- function(train_data, test_data, encoding_dim = 64, epochs = 50) {
  # Define the input shape
  input_shape <- dim(train_data)[2]
  
  # Define the encoder network
  encoder <- keras_model_sequential() %>%
    layer_dense(units = 128, activation = "relu", input_shape = input_shape) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = as.integer(encoding_dim), activation = "relu") 
  
  # Define the decoder network
  decoder <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = "relu", input_shape = encoding_dim) %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = input_shape, activation = "sigmoid")
  
  # Define the autoencoder model
  autoencoder <- keras_model_sequential() %>%
    add(encoder) %>%
    add(decoder)
  
  # Compile the model
  autoencoder %>% compile(
    optimizer = "adam",
    loss = "mse"
  )
  
  # Train the model
  autoencoder %>% fit(
    train_data, train_data,
    epochs = epochs,
    batch_size = 32,
    shuffle = TRUE,
    validation_data = list(test_data, test_data)
  )
  
  # Extract the encoded data
  encoded_data <- encoder %>% predict(train_data)
  
  # Extract the decoded data
  decoded_data <- decoder %>% predict(encoded_data)
  
  # Return the encoded and decoded data as a list
  return(list(encoded_data = encoded_data, decoded_data = decoded_data))
}

# Load necessary packages
library(keras)

# Generate some example SNP data
train_snp <- matrix(sample(c(0, 1, 2), size = 1000*100, replace = TRUE), nrow = 1000)

# Split data into training and test sets
train_data <- train_snp[1:800, ]
test_data <- train_snp[801:1000, ]

# Run the autoencoder function
autoencoder_result <- autoencoder_func(train_data, test_data, encoding_dim = 10, epochs = 50)

# Extract the encoded and decoded data
encoded_data <- autoencoder_result$encoded_data
decoded_data <- autoencoder_result$decoded_data



# sol_VL
autoencoder_sol_VL <- autoencoder_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- autoencoder_sol_VL$val_predicted
train_predicted <- autoencoder_sol_VL$predicted_marker
x <- (autoencoder_sol_VL$predicted_marker)
summary(sol_VL_test_phenotype)
summary(test_predicted)

cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

sol_VL_train_test_autoencoder <- cbind(sol_VL_test_phenotype, test_predicted)
