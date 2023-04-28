library(ggplot2)
# Filtering out the significant markers from GWAS
sol_Mo_gwas_markers <- unlist(as.vector(sol_Mo_gwas[,1]))
hist(sol_Mo$sol_Mo, xlab = "Probability of finding moderate P", main = "sol_Mo")
# Filtering out the significant markers from GWAS
sol_Mo_gwas_markers <- unlist(as.vector(sol_Mo_gwas[,1]))

##### ASSUMPTION - the missing makers are homozygous
# Subsetting markers - replacing -9 with 0
sol_Mo_markers <- SNP_markers[sol_Mo_gwas_markers]
sol_Mo_markers <- replace(sol_Mo_markers, sol_Mo_markers == -9,0)

# Set the number of iterations to run the loop
num_iter <- 30

# Correlation
train_corr_sol_Mo <- vector("numeric", num_iter)
test_corr_sol_Mo <- vector("numeric", num_iter)


# Create empty vectors to store the predicted values for each iteration
train_preds_sol_Mo <- matrix(NA, nrow = nrow(sol_Mo), ncol = num_iter)
test_preds_sol_Mo <- matrix(NA, nrow = nrow(sol_Mo), ncol = num_iter)

# Initialize lists to store test predictions
test_preds_sol_Mo_list <- vector("list", num_iter)

# Initialize empty matrices to store the test predictions for each phenotype
sol_Mo_pred_matrix <- matrix(NA, nrow = nrow(sol_Mo), 1)
i = 1
# Run the loop
for (i in 1:num_iter) {
  print(i)
  error_flag <- TRUE
  while (error_flag) {
    # Randomly sample markers and phenotypes for the training set
    sol_Mo_train <- as.matrix(sample(1:nrow(sol_Mo_markers), 1035))
    sol_Mo_test <- setdiff(1:nrow(sol_Mo_markers), sol_Mo_train)
    sol_Mo_train_phenotype <- as.matrix(sol_Mo[sol_Mo_train, ])
    sol_Mo_train_marker <- as.matrix(sol_Mo_markers[sol_Mo_train, ], K = NULL)
    # Removing NAs from training and putting them into testing
    # nas <- which(is.na(sol_Mo_train_phenotype))
    # sol_Mo_train_phenotype <- as.matrix(sol_Mo_train_phenotype[-nas, ])
    # sol_Mo_train_marker <- as.matrix(sol_Mo_train_marker[-nas, ])
    # sol_Mo_train <- sol_Mo_train[-nas,]
    # sol_Mo_test <- sort(c(sol_Mo_test, nas))
    sol_Mo_test_phenotype <- as.matrix(sol_Mo[sol_Mo_test, ])
    sol_Mo_test_marker <- as.matrix(sol_Mo_markers[sol_Mo_test, ], K = NULL)
    
    
    
    # Wrap the DNN_func() call inside a tryCatch() block
    error_flag <- FALSE
    tryCatch({
      DNN_sol_Mo <- DNN_func(sol_Mo_train_phenotype, sol_Mo_train_marker, SNP_markers)
      
      # Extract the predicted values for the training and test sets
      train_preds_sol_Mo[1:nrow(DNN_sol_Mo$train_predicted), i] <- DNN_sol_Mo$train_predicted
      test_preds_sol_Mo[1:nrow(DNN_sol_Mo$val_predicted), i] <- DNN_sol_Mo$val_predicted
      
      #Calculate correlation without NA values
      train_corr_sol_Mo[i] <- cor(sol_Mo_train_phenotype, train_preds_sol_Mo[1:nrow(DNN_sol_Mo$train_predicted), i], use = "pairwise.complete.obs")
      test_corr_sol_Mo[i] <- cor(sol_Mo_test_phenotype, test_preds_sol_Mo[1:nrow(DNN_sol_Mo$val_predicted), i], use = "pairwise.complete.obs")
      
      
      
      # Fill the matrices with the test predictions based on the test indices
      sol_Mo_pred_matrix[sol_Mo_test, 1] <- DNN_sol_Mo$val_predicted
      
      # Store the filled matrices in the respective lists
      test_preds_sol_Mo_list[[i]] <- sol_Mo_pred_matrix
    }, error = function(e) {
      error_flag <- TRUE
    })
    
    # # If an error occurred during DNN_func(), skip the rest of the current iteration
    # if (error_flag) {
    #   next
    # }
  }

}

# Calculate mean correlation
#mean_train_corr_sol_Mo <- mean(train_corr_sol_Mo)
#mean_test_corr_sol_Mo <- mean(test_corr_sol_Mo)

# saveRDS(train_corr_sol_Mo,"train_corr_sol_Mo.RDS")
# saveRDS(test_corr_sol_Mo,"test_corr_sol_Mo.RDS")
# saveRDS(test_preds_sol_Mo_list,"test_preds_sol_Mo_list.RDS")
# saveRDS(DNN_sol_Mo, "DNN_sol_Mo.RDS")

# train_corr_sol_Mo <- readRDS("train_corr_sol_Mo.RDS")
# test_corr_sol_Mo <- readRDS("test_corr_sol_Mo.RDS")
# test_preds_sol_Mo_list <- readRDS("test_preds_sol_Mo_list.RDS")
# DNN_sol_Mo <- readRDS("DNN_sol_Mo.RDS")

# Combine the list of matrices into a single data frame for each phenotype
test_preds_sol_Mo_df <- do.call(cbind, test_preds_sol_Mo_list)

# Set the column names for each data frame
colnames(test_preds_sol_Mo_df) <- paste0("Iteration_", 1:num_iter)
rownames(test_preds_sol_Mo_df) <- rownames(sol_Mo)

########## Correlation
# Calculate mean and standard deviation for train and test correlations
train_mean <- mean(train_corr_sol_Mo, na.rm = TRUE)
train_sd <- sd(train_corr_sol_Mo, na.rm = TRUE)
test_mean <- mean(test_corr_sol_Mo, na.rm = TRUE)
test_sd <- sd(test_corr_sol_Mo, na.rm = TRUE)

# Create a dataframe with the summary statistics
summary_data <- data.frame(
  Category = c("Train", "Test"),
  Mean = c(train_mean, test_mean),
  StdDev = c(train_sd, test_sd)
)

# Create a bar plot with error bars
w <- ggplot(summary_data, aes(x = Category, y = Mean, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
  geom_errorbar(aes(ymin = Mean - StdDev, ymax = Mean + StdDev), 
                position = position_dodge(0.6), width = 0.2) +
  geom_text(aes(label = round(Mean, 4)), 
            position = position_dodge(width = 0.6), vjust = -0.5) +
  labs(title = "Mean and Standard Deviation of Train and Test Correlations",
       x = "Category",
       y = "Mean Correlation") +
  theme_minimal() +
  theme(legend.position = "none")




############# Best and Worst Accessions for SD
# Remove NAs from the data frame
#test_preds_sol_Mo_df <- as.data.frame(rowMeans(test_preds_sol_Mo_df, na.rm =T))
#test_preds_sol_Mo_df <- test_preds_sol_Mo_df[complete.cases(test_preds_sol_Mo_df),]

# Count the number of non-NA values in each row
non_na_counts <- length(apply(!is.na(test_preds_sol_Mo_df), 1, sum))
non_na_counts

# Calculate the mean and standard error for each row, ignoring NA values
means <- apply(test_preds_sol_Mo_df, 1, function(x) mean(x, na.rm = TRUE))
std_errors <- apply(test_preds_sol_Mo_df, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Create a data frame with the calculated values
data <- data.frame(Sample = 1:nrow(test_preds_sol_Mo_df), Mean = means, StdError = std_errors)

# Sort the data frame by the standard error
data <- data[order(data$StdError),]

# Select the best 25 samples with the lowest standard error and the worst 25 samples with the highest standard error
selected_data <- rbind(head(data, 25), tail(data, 25))

# Add a column to identify best and worst samples
selected_data$Type <- c(rep("Best", 25), rep("Worst", 25))

# Create the bar plot with error bars
x <- ggplot(selected_data, aes(x = reorder(Sample, -StdError), y = Mean, fill = Type)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Mean - StdError, ymax = Mean + StdError), width = 0.2) +
  xlab("Sample") + ylab("Mean") + ggtitle("Best 25 and Worst 25 Samples by Standard Error") +
  scale_x_discrete(labels = selected_data$Sample) +
  scale_fill_manual(values = c("Best" = "steelblue", "Worst" = "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))





############# Best and Worst Accessions for Phenotype

# Create a data frame from the original sol_Mo data
sol_Mo_df <- data.frame(Sample = rownames(sol_Mo), sol_Mo = sol_Mo[, "sol_Mo"])
# Merge the predicted data with the original data
merged_data <- merge(data, sol_Mo_df, by = "Sample")
# Calculate the absolute difference between the predicted and original values
merged_data$Difference <- abs(merged_data$Mean - merged_data$sol_Mo)
# Remove rows with NA values in the Difference column
merged_data <- merged_data[!is.na(merged_data$Difference),]
# Sort the samples by the absolute difference
sorted_data <- merged_data[order(merged_data$Difference), ]
# Select the top 25 and worst 25 samples
top_25 <- sorted_data[1:25, ]
worst_25 <- sorted_data[(nrow(sorted_data) - 24):nrow(sorted_data), ]
# Combine the top 25 and worst 25 samples
selected_data <- rbind(top_25, worst_25)
# Add a column to indicate whether a sample is among the best or worst
selected_data$Category <- c(rep("Best", 25), rep("Worst", 25))
selected_data$Sample <- as.character(selected_data$Sample)
# Filter the best and worst samples
best_samples <- selected_data[selected_data$Category == "Best", ]
worst_samples <- selected_data[selected_data$Category == "Worst", ]

# Create a bar plot for the best samples
y <- ggplot(best_samples, aes(x = Sample, y = Difference)) +
  geom_bar(stat = "identity", fill = "dodgerblue") +
  labs(title = "Best Samples: Difference between Predicted and Original Values",
       x = "Sample",
       y = "Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Create a bar plot for the worst samples
z <- ggplot(worst_samples, aes(x = Sample, y = Difference)) +
  geom_bar(stat = "identity", fill = "red") +
  labs(title = "Worst Samples: Difference between Predicted and Original Values",
       x = "Sample",
       y = "Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

###### Arranging the plots
library(gridExtra)

arranged_plot <- grid.arrange(w, x,
             y, z,
             nrow = 2, ncol = 2)

# Save the arranged plot as a PNG file
ggsave("sol_Mo_DNN_no9to0.png", arranged_plot, width = 15, height = 10, units = "in", dpi = 300)


table(is.na(means))
