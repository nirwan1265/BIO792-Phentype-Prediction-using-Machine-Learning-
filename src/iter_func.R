# Loading the phenotype data
# Change working directory
dir_pheno <- "~/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/phenotype_data/"

# Phenotypes
phenotypes <- read.csv(paste0(dir_pheno,"Sorghum_allphospho_africa.csv"))

# Total phosphorus
stp <- phenotypes["stp10"]

# Solubility
PBR <- phenotypes["PBR1"]
#hist(lab[,1])

# Phosphorus retention
PNZ <- phenotypes["PNZ1"]
#hist(sol_Mo[,1])

# Alkaline soil P
sol_Mo <- phenotypes["sol_Mo"]


####### Filtering process
### SNP selection using GWAS MLM 
# Loading GWAS files and only selecting the significant snps
dir_gwas <- "~/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/gwas/"
stp_gwas <- vroom(paste0(dir_gwas,"stp10.txt")) %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
PBR_gwas <- vroom(paste0(dir_gwas,"PBR1.txt")) %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
PNZ_gwas <- vroom(paste0(dir_gwas,"PNZ1.txt")) %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
sol_Mo_gwas <- vroom(paste0(dir_gwas,"sol_Mo.txt")) %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)

####### Loading the genotype file (MAF)
# Imputed using- 1.) nearest neighbor searches using Viterbi algorithm and then LDKiNN
dir_geno <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/v3/imputed/"
SNP_markers <- vroom(paste0(dir_geno,"allchrom_MAF_sorghum.txt"))
# Sanity check
#SNP_markers[1:10,1:10]

# Filtering out the significant markers from GWAS
stp_gwas_markers <- unlist(as.vector(stp_gwas[,1]))
PBR_gwas_markers <- unlist(as.vector(PBR_gwas[,1]))
PNZ_gwas_markers <- unlist(as.vector(PNZ_gwas[,1]))
sol_Mo_gwas_markers <- unlist(as.vector(sol_Mo_gwas[,1]))

##### ASSUMPTION - the missing makers are homozygous
# Subsetting markers - replacing -9 with 0
stp_markers <- SNP_markers[stp_gwas_markers]
stp_markers <- replace(stp_markers, stp_markers == -9,0)

PBR_markers <- SNP_markers[PBR_gwas_markers]
PBR_markers <- replace(PBR_markers, PBR_markers == -9,0)

PNZ_markers <- SNP_markers[PNZ_gwas_markers]
PNZ_markers <- replace(PNZ_markers, PNZ_markers == -9,0)

sol_Mo_markers <- SNP_markers[sol_Mo_gwas_markers]
sol_Mo_markers <- replace(sol_Mo_markers, sol_Mo_markers == -9,0)

# Set the number of iterations to run the loop
num_iter <- 100

# Correlation
train_corr_stp <- vector("numeric", num_iter)
test_corr_stp <- vector("numeric", num_iter)

train_corr_PBR <- vector("numeric", num_iter)
test_corr_PBR <- vector("numeric", num_iter)

train_corr_PNZ <- vector("numeric", num_iter)
test_corr_PNZ <- vector("numeric", num_iter)

train_corr_sol_Mo <- vector("numeric", num_iter)
test_corr_sol_Mo <- vector("numeric", num_iter)

# Create empty vectors to store the predicted values for each iteration
train_preds_stp <- matrix(NA, nrow = nrow(stp), ncol = num_iter)
test_preds_stp <- matrix(NA, nrow = nrow(stp), ncol = num_iter)

train_preds_PBR <- matrix(NA, nrow = nrow(PBR), ncol = num_iter)
test_preds_PBR <- matrix(NA, nrow = nrow(PBR), ncol = num_iter)

train_preds_PNZ <- matrix(NA, nrow = nrow(PNZ), ncol = num_iter)
test_preds_PNZ <- matrix(NA, nrow = nrow(PNZ), ncol = num_iter)

train_preds_sol_Mo <- matrix(NA, nrow = nrow(sol_Mo), ncol = num_iter)
test_preds_sol_Mo <- matrix(NA, nrow = nrow(sol_Mo), ncol = num_iter)

# Initialize lists to store test predictions
test_preds_stp_list <- vector("list", num_iter)
test_preds_PBR_list <- vector("list", num_iter)
test_preds_PNZ_list <- vector("list", num_iter)
test_preds_sol_Mo_list <- vector("list", num_iter)


# Run the loop
for (i in 1:num_iter) {
  
  # Randomly sample markers and phenotypes for the training set
  stp_train <- as.matrix(sample(1:nrow(stp_markers), 1035))
  stp_test <- setdiff(1:nrow(stp_markers), stp_train)
  stp_train_phenotype <- as.matrix(stp[stp_train, ])
  
  stp_train_marker <- as.matrix(stp_markers[stp_train, ], K = NULL)
  #Removing NAs from training and putting them into testing
  nas <- which(is.na(stp_train_phenotype))
  stp_train_phenotype <- as.matrix(stp_train_phenotype[-nas, ])
  stp_train_marker <- as.matrix(stp_train_marker[-nas, ])
  stp_train <- stp_train[-nas,]
  stp_test <- sort(c(stp_test, nas))
  stp_test_phenotype <- as.matrix(stp[stp_test, ])
  stp_test_marker <- as.matrix(stp_markers[stp_test, ], K = NULL)
  print(nrow(stp_train_phenotype))
  # Randomly sample markers and phenotypes for the training set
  # PBR_train <- as.matrix(sample(1:nrow(PBR_markers), 1035))
  # PBR_test <- setdiff(1:nrow(PBR_markers), PBR_train)
  # PBR_train_phenotype <- as.matrix(PBR[PBR_train, ])
  # PBR_train_marker <- as.matrix(PBR_markers[PBR_train, ], K = NULL)
  # # Removing NAs from training and putting them into testing
  # nas <- which(is.na(PBR_train_phenotype))
  # PBR_train_phenotype <- as.matrix(PBR_train_phenotype[-nas, ])
  # PBR_train_marker <- as.matrix(PBR_train_marker[-nas, ])
  # PBR_train <- PBR_train[-nas,]
  # PBR_test <- sort(c(PBR_test, nas))
  # PBR_test_phenotype <- as.matrix(PBR[PBR_test, ])
  # PBR_test_marker <- as.matrix(PBR_markers[PBR_test, ], K = NULL)
  # 
  # # Randomly sample markers and phenotypes for the training set
  # PNZ_train <- as.matrix(sample(1:nrow(PNZ_markers), 1035))
  # PNZ_test <- setdiff(1:nrow(PNZ_markers), PNZ_train)
  # PNZ_train_phenotype <- as.matrix(PNZ[PNZ_train, ])
  # PNZ_train_marker <- as.matrix(PNZ_markers[PNZ_train, ], K = NULL)
  # # Removing NAs from training and putting them into testing
  # nas <- which(is.na(PNZ_train_phenotype))
  # PNZ_train_phenotype <- as.matrix(PNZ_train_phenotype[-nas, ])
  # PNZ_train_marker <- as.matrix(PNZ_train_marker[-nas, ])
  # PNZ_train <- PNZ_train[-nas,]
  # PNZ_test <- sort(c(PNZ_test, nas))
  # PNZ_test_phenotype <- as.matrix(PNZ[PNZ_test, ])
  # PNZ_test_marker <- as.matrix(PNZ_markers[PNZ_test, ], K = NULL)
  # 
  # # Randomly sample markers and phenotypes for the training set
  # sol_Mo_train <- as.matrix(sample(1:nrow(sol_Mo_markers), 1035))
  # sol_Mo_test <- setdiff(1:nrow(sol_Mo_markers), sol_Mo_train)
  # sol_Mo_train_phenotype <- as.matrix(sol_Mo[sol_Mo_train, ])
  # sol_Mo_train_marker <- as.matrix(sol_Mo_markers[sol_Mo_train, ], K = NULL)
  # sol_Mo_test_phenotype <- as.matrix(sol_Mo[sol_Mo_test, ])
  # sol_Mo_test_marker <- as.matrix(sol_Mo_markers[sol_Mo_test, ], K = NULL)
  
  # Fit the SVM model to the training data
  SVM_stp <- SVM_func(stp_train_phenotype, stp_train_marker, stp_test_marker)
  # SVM_PBR <- SVM_func(PBR_train_phenotype, PBR_train_marker, PBR_test_marker)
  # SVM_PNZ <- SVM_func(PNZ_train_phenotype, PNZ_train_marker, PNZ_test_marker)
  # SVM_sol_Mo <- SVM_func(sol_Mo_train_phenotype, sol_Mo_train_marker, sol_Mo_test_marker)
  # 
  # Extract the predicted values for the training and test sets
  train_preds_stp[1:nrow(SVM_stp$train_predicted), i] <- SVM_stp$train_predicted
  test_preds_stp[1:nrow(SVM_stp$val_predicted), i] <- SVM_stp$val_predicted
  # train_preds_PBR[1:nrow(SVM_PBR$train_predicted), i] <- SVM_PBR$train_predicted
  # test_preds_PBR[1:nrow(SVM_PBR$val_predicted), i] <- SVM_PBR$val_predicted
  # train_preds_PNZ[1:nrow(SVM_PNZ$train_predicted), i] <- SVM_PNZ$train_predicted
  # test_preds_PNZ[1:nrow(SVM_PNZ$val_predicted), i] <- SVM_PNZ$val_predicted
  # train_preds_sol_Mo[1:nrow(SVM_sol_Mo$train_predicted), i] <- SVM_sol_Mo$train_predicted
  # test_preds_sol_Mo[1:nrow(SVM_sol_Mo$val_predicted), i] <- SVM_sol_Mo$val_predicted
  # 
  
  # Calculate correlation without NA values
  train_corr_stp[i] <- cor(stp_train_phenotype, train_preds_stp[1:nrow(SVM_stp$train_predicted), i], use = "pairwise.complete.obs")
  test_corr_stp[i] <- cor(stp_test_phenotype, test_preds_stp[1:nrow(SVM_stp$val_predicted), i], use = "pairwise.complete.obs")
  # train_corr_PBR[i] <- cor(PBR_train_phenotype, train_preds_PBR[1:nrow(SVM_PBR$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_PBR[i] <- cor(PBR_test_phenotype, test_preds_PBR[1:nrow(SVM_PBR$val_predicted), i], use = "pairwise.complete.obs")
  # train_corr_PNZ[i] <- cor(PNZ_train_phenotype, train_preds_PNZ[1:nrow(SVM_PNZ$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_PNZ[i] <- cor(PNZ_test_phenotype, test_preds_PNZ[1:nrow(SVM_PNZ$val_predicted), i], use = "pairwise.complete.obs")
  # train_corr_sol_Mo[i] <- cor(sol_Mo_train_phenotype, train_preds_sol_Mo[1:nrow(SVM_sol_Mo$train_predicted), i], use = "pairwise.complete.obs")
  # test_corr_sol_Mo[i] <- cor(sol_Mo_test_phenotype, test_preds_sol_Mo[1:nrow(SVM_sol_Mo$val_predicted), i], use = "pairwise.complete.obs")
  # 
  # Initialize empty matrices to store the test predictions for each phenotype
  stp_pred_matrix <- matrix(NA, nrow = nrow(stp), 1)
  # PBR_pred_matrix <- matrix(NA, nrow = nrow(PBR), 1)
  # PNZ_pred_matrix <- matrix(NA, nrow = nrow(PNZ), 1)
  # sol_Mo_pred_matrix <- matrix(NA, nrow = nrow(sol_Mo), 1)
  # 
  # Fill the matrices with the test predictions based on the test indices
  stp_pred_matrix[stp_test, 1] <- SVM_stp$val_predicted
  # PBR_pred_matrix[PBR_test, 1] <- SVM_PBR$val_predicted
  # PNZ_pred_matrix[PNZ_test, 1] <- SVM_PNZ$val_predicted
  # sol_Mo_pred_matrix[sol_Mo_test, 1] <- SVM_sol_Mo$val_predicted
  # 
  # Store the filled matrices in the respective lists
  test_preds_stp_list[[i]] <- stp_pred_matrix
  # test_preds_PBR_list[[i]] <- PBR_pred_matrix
  # test_preds_PNZ_list[[i]] <- PNZ_pred_matrix
  # test_preds_sol_Mo_list[[i]] <- sol_Mo_pred_matrix
}


# Calculate mean correlation
mean_train_corr_stp <- mean(train_corr_stp)
mean_test_corr_stp <- mean(test_corr_stp)
# mean_train_corr_PBR <- mean(train_corr_PBR)
# mean_test_corr_PBR <- mean(test_corr_PBR)
# mean_train_corr_PNZ <- mean(train_corr_PNZ)
# mean_test_corr_PNZ <- mean(test_corr_PNZ)
# mean_train_corr_sol_Mo <- mean(train_corr_sol_Mo)
# mean_test_corr_sol_Mo <- mean(test_corr_sol_Mo)

# Combine the list of matrices into a single data frame for each phenotype
test_preds_stp_df <- do.call(cbind, test_preds_stp_list)
# test_preds_PBR_df <- do.call(cbind, test_preds_PBR_list)
# test_preds_PNZ_df <- do.call(cbind, test_preds_PNZ_list)
# test_preds_sol_Mo_df <- do.call(cbind, test_preds_sol_Mo_list)

# Set the column names for each data frame
colnames(test_preds_stp_df) <- paste0("Iteration_", 1:37)
# colnames(test_preds_PBR_df) <- paste0("Iteration_", 1:num_iter)
# colnames(test_preds_PNZ_df) <- paste0("Iteration_", 1:num_iter)
# colnames(test_preds_sol_Mo_df) <- paste0("Iteration_", 1:num_iter)

# Calculate the mean phenotype per row without NAs for each data frame
mean_test_preds_stp <- rowMeans(test_preds_stp_df, na.rm = TRUE)
# mean_test_preds_PBR <- rowMeans(test_preds_PBR_df, na.rm = TRUE)
# mean_test_preds_PNZ <- rowMeans(test_preds_PNZ_df, na.rm = TRUE)
# mean_test_preds_sol_Mo <- rowMeans(test_preds_sol_Mo_df, na.rm = TRUE)

# Create data frames with just the mean values for each phenotype
mean_test_preds_stp_df <- data.frame(Mean_Prediction = mean_test_preds_stp)
# mean_test_preds_PBR_df <- data.frame(Mean_Prediction = mean_test_preds_PBR)
# mean_test_preds_PNZ_df <- data.frame(Mean_Prediction = mean_test_preds_PNZ)
# mean_test_preds_sol_Mo_df <- data.frame(Mean_Prediction = mean_test_preds_sol_Mo)




a <- test_preds_stp_df
b <- mean_test_preds_stp_df
c <- train_corr_stp
d <- test_corr_stp



x <- test_preds_stp_df
y <- mean_test_preds_stp_df
trc <- train_corr_stp
tc <- test_corr_stp



rm(list = setdiff(ls(), "SNP_markers"))
