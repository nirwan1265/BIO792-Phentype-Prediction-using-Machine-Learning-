# Remove the last two layers from the DNN model
DNN_results <- DNN_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)
DNN_results <- DNN_sol_VL
DNN_feature_extractor <- keras_model(inputs = DNN_results$model$input, outputs = DNN_results$model$layers[[length(DNN_results$model$layers) - 2]]$output)


# Extract features from the training and testing genotypes
train_geno_features <- predict(DNN_feature_extractor, sol_VL_train_marker)
test_geno_features <- predict(DNN_feature_extractor, sol_VL_test_marker)

#Train the RNN model with the extracted feature
RNN_results <- RNN_func(sol_VL_train_phenotype, train_geno_features, test_geno_features)


test_predicted <- RNN_results$val_pred

train_predicted <- RNN_results$train_pred


cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

summary(sol_VL_test_phenotype)
summary(test_predicted)
summary(sol_VL_test_phenotype)

plot(sol_VL_test_phenotype,test_predicted,xlab = "True Values", ylab = "Predicted Values", main = "True vs Predicted Values")


sol_VL_train_test_rnn <- cbind(sol_VL_test_phenotype, test_predicted)
colnames(sol_VL_train_test_rnn) <- c("test","train")
