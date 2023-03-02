library(vroom)
# Change working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data")
# Phenotypes
phenotypes <- read.csv("Sorghum_allphospho_africa.csv")
# Total phosphorus
tot <- phenotypes["tot"]
hist(tot[,1])

# Solubility
lab <- phenotypes["lab"]
hist(lab[,1])

# Phosphorus retention
sol_VL <- phenotypes["sol_VL"]
hist(sol_VL[,1])


# Loading GWAS files and only selecting the significant snps
tot_gwas <- vroom("tot_LMM.txt") %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
lab_gwas <- vroom("lab_LMM.txt") %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
sol_VL_gwas <- vroom("sol_VL_LMM.txt") %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
#tot_gwas_trial <- vroom("tot_LMM.txt") %>% select(rs,p_wald)

# Loading the genotype file (MAF)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/v3/imputed")
#Imputed using LDKiNN
SNP_markers <- vroom("combined.txt")
#SNP_markers[1:100,1:100]
#SNP_markers <- vroom("allchrom_africa_filtered.MAF.txt")
#saveRDS(SNP_markers,"SNP_markers.RDS")
#SNP_markers <- readRDS("SNP_markers.RDS")
SNP_markers[1:10,1:10]

# Getting the significant markers
tot_gwas_markers <- unlist(as.vector(tot_gwas[,1]))
lab_gwas_markers <- unlist(as.vector(lab_gwas[,1]))
sol_VL_gwas_markers <- unlist(as.vector(sol_VL_gwas[,1]))


# Subsetting markers - replacing -9 with 0
tot_markers <- SNP_markers[tot_gwas_markers]
tot_markers <- replace(tot_markers, tot_markers == -9,0)

lab_markers <- SNP_markers[lab_gwas_markers]
lab_markers <- replace(lab_markers, lab_markers == -9,0)

sol_VL_markers <- SNP_markers[sol_VL_gwas_markers]
sol_VL_markers <- replace(sol_VL_markers, sol_VL_markers == -9,0)

# Subsetting Testing and Training data for markers
set.seed(123)

tot_train <- as.matrix(sample(1:nrow(tot_markers), 1035))
tot_test <- setdiff(1:nrow(tot_markers), tot_train)

lab_train <- as.matrix(sample(1:nrow(lab_markers), 1035))
lab_test <- setdiff(1:nrow(lab_markers), lab_train)

sol_VL_train <- as.matrix(sample(1:nrow(sol_VL_markers), 1035))
sol_VL_test <- setdiff(1:nrow(sol_VL_markers), sol_VL_train)


# Subsetting Testing and Training data for phenotypes and markers

tot_train_phenotype <- as.matrix(tot[tot_train, ])
tot_train_marker <- as.matrix(tot_markers[tot_train, ], K = NULL)
# Removing NAs from training and putting them into testing
nas <- which(is.na(tot_train_phenotype))
tot_train_phenotype <- as.matrix(tot_train_phenotype[-nas, ])
tot_train_marker <- as.matrix(tot_train_marker[-nas, ])
tot_train <- tot_train[-nas,]
tot_test <- sort(c(tot_test, nas))
tot_test_phenotype <- as.matrix(tot[tot_test, ])
tot_test_marker <- as.matrix(tot_markers[tot_test, ], K = NULL)


# # Subsetting Testing and Training data for phenotypes and markers
# tot_train_phenotype <- as.matrix(tot[tot_train, ])
# # removing NA and updating them to the test data
# na_rows <- which(is.na(tot_train_phenotype))
# tot_test <- sort(c(tot_test,na_rows))
# tot_train <- setdiff(1:nrow(tot_markers), tot_test)
# tot_train_phenotype <- as.matrix(tot_train_phenotype[!is.na(tot_train_phenotype)])
# tot_train_marker <- as.matrix(tot_markers[tot_train, ], K = NULL)
# tot_test_phenotype <- as.matrix(tot[tot_test, ])
# tot_test_marker <- as.matrix(tot_markers[tot_test, ], K = NULL)

lab_train_phenotype <- as.matrix(lab[lab_train, ])
lab_train_marker <- as.matrix(lab_markers[lab_train, ], K = NULL)
# Removing NAs from training and putting them into testing
nas <- which(is.na(lab_train_phenotype))
lab_train_phenotype <- as.matrix(lab_train_phenotype[-nas, ])
lab_train_marker <- as.matrix(lab_train_marker[-nas, ])
lab_train <- lab_train[-nas,]
lab_test <- sort(c(lab_test, nas))
lab_test_phenotype <- as.matrix(lab[lab_test, ])
lab_test_marker <- as.matrix(lab_markers[lab_test, ], K = NULL)


sol_VL_train_phenotype <- as.matrix(sol_VL[sol_VL_train, ])
sol_VL_train_marker <- as.matrix(sol_VL_markers[sol_VL_train, ], K = NULL)
# Removing NAs from training and putting them into testing
nas <- which(is.na(sol_VL_train_phenotype))

sol_VL_train_phenotype <- as.matrix(sol_VL_train_phenotype[-nas, ])
sol_VL_train_marker <- as.matrix(sol_VL_train_marker[-nas, ])
sol_VL_train <- sol_VL_train[-nas,]
sol_VL_test <- sort(c(sol_VL_test, nas))
sol_VL_test_phenotype <- as.matrix(sol_VL[sol_VL_test, ])
sol_VL_test_marker <- as.matrix(sol_VL_markers[sol_VL_test, ], K = NULL)



# sol_VL_train_phenotype <- as.matrix(sol_VL[sol_VL_train, ])
# # removing NA and updating them to the test data
# na_rows <- which(is.na(sol_VL_train_phenotype))
# sol_VL_test <- sort(c(sol_VL_test,na_rows))
# sol_VL_train <- setdiff(1:nrow(sol_VL_markers), sol_VL_test)
# sol_VL_train_phenotype <- as.matrix(sol_VL_train_phenotype[!is.na(sol_VL_train_phenotype)])
# sol_VL_train_marker <- as.matrix(sol_VL_markers[sol_VL_train, ], K = NULL)
# sol_VL_test_phenotype <- as.matrix(sol_VL[sol_VL_test, ])
# sol_VL_test_marker <- as.matrix(sol_VL_markers[sol_VL_test, ], K = NULL)


# Training the model
model_tot <- mixed.solve(y = tot_train_phenotype, Z = tot_train_marker)

model_lab <- mixed.solve(y = lab_train_phenotype, Z = lab_train_marker)

model_sol_VL <- mixed.solve(y = sol_VL_train_phenotype, Z = sol_VL_train_marker)


# Marker effects
marker_effects_tot <- as.matrix(model_tot$u)

marker_effects_lab <- as.matrix(model_lab$u)

marker_effects_sol_VL <- as.matrix(model_sol_VL$u)


# Running BLUE
BLUE_tot <- as.vector(model_tot$beta)

BLUE_lab <- as.vector(model_lab$beta)

BLUE_sol_VL <- as.vector(model_sol_VL$beta)


# Prediction
predict_train_tot <- as.matrix(tot_train_marker) %*% marker_effects_tot
predict_test_tot <- as.matrix(tot_test_marker) %*% marker_effects_tot
predict_train_result_tot <- as.vector((predict_train_tot[,1])+BLUE_tot)
predict_test_result_tot <- as.vector((predict_test_tot[,1])+BLUE_tot)

predict_train_lab <- as.matrix(lab_train_marker) %*% marker_effects_lab
predict_test_lab <- as.matrix(lab_test_marker) %*% marker_effects_lab
predict_train_result_lab <- as.vector((predict_train_lab[,1])+BLUE_lab)
predict_test_result_lab <- as.vector((predict_test_lab[,1])+BLUE_lab)

predict_train_sol_VL <- as.matrix(sol_VL_train_marker) %*% marker_effects_sol_VL
predict_test_sol_VL <- as.matrix(sol_VL_test_marker) %*% marker_effects_sol_VL
predict_train_result_sol_VL <- as.vector((predict_train_sol_VL[,1])+BLUE_sol_VL)
predict_test_result_sol_VL <- as.vector((predict_test_sol_VL[,1])+BLUE_sol_VL)


# Correlation
corr_tot_train <- cor.test(as.vector(tot_train_phenotype), predict_train_result_tot, use = "spearman")
corr_tot_test <- cor(as.vector(tot_test_phenotype), predict_test_result_tot, use = "complete")

corr_lab_train <- cor.test(as.vector(lab_train_phenotype), predict_train_result_lab, use = "spearman")
corr_lab_test <- cor(as.vector(lab_test_phenotype), predict_test_result_lab, use = "complete")

corr_sol_VL_train <- cor.test(as.vector(sol_VL_train_phenotype), predict_train_result_sol_VL, use = "spearman")
corr_sol_VL_test <- cor(as.vector(sol_VL_test_phenotype), predict_test_result_sol_VL, use = "complete")


# Visualize
plot(tot_test_phenotype, predict_test_result_tot)
abline(lm(tot_test_phenotype~predict_test_result_tot), col = "red")

plot(lab_test_phenotype, predict_test_result_lab)
plot(sol_VL_test_phenotype, predict_test_result_sol_VL)
abline(lm(sol_VL_test_phenotype~ predict_test_result_sol_VL), col = "red")

# Visualizing side by side
tot_train_test <- cbind(tot_test_phenotype, predict_test_result_tot)
lab_train_test <- cbind(lab_test_phenotype, predict_test_result_lab)
sol_VL_train_test <- cbind(sol_VL_test_phenotype, predict_test_result_sol_VL)

x <- cbind(sol_VL_train_test,sol_VL_train_test_DNN)
x <- x[,-3]
colnames(x) <- c("soilP","rrBLUP","CNN")
x <- as.data.frame(x)
str(x)
x$CNN <- x$CNN * 2





# SNP interaction data 
snp_interactions is the matrix of interactions
A <- A.mat(snp_interactions, model = "Additive")
result <- mixed.solve(y = y, X = X, K = A)

