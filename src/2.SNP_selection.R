# Loading the phenotype data
# Change working directory
dir_pheno <- "~/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/phenotype_data/"

# Phenotypes
phenotypes <- read.csv(paste0(dir_pheno,"Sorghum_allphospho_africa.csv"))

# Total phosphorus
stp <- phenotypes["stp10"]

# Solubility
sol_Mo <- phenotypes["sol_Mo"]
#hist(lab[,1])

# Phosphorus retention
PNZ <- phenotypes["PNZ1"]
#hist(sol_Mo[,1])

# Alkaline soil P
POL <- phenotypes["POL1"]


####### Filtering process
### SNP selection using GWAS MLM 
# Loading GWAS files and only selecting the significant snps
dir_gwas <- "~/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/gwas/"
stp_gwas <- vroom(paste0(dir_gwas,"stp10.txt")) %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
sol_Mo_gwas <- vroom(paste0(dir_gwas,"sol_Mo.txt")) %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
PNZ_gwas <- vroom(paste0(dir_gwas,"PNZ1.txt")) %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)
POL_gwas <- vroom(paste0(dir_gwas,"POL1.txt")) %>% select(rs,p_wald) %>% filter(p_wald <= 0.05)

####### Loading the genotype file (MAF)
# Imputed using- 1.) nearest neighbor searches using Viterbi algorithm and then LDKiNN
dir_geno <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/v3/imputed/"
SNP_markers <- vroom(paste0(dir_geno,"allchrom_MAF_sorghum.txt"))
# Sanity check
SNP_markers[1:10,1:10]

# Filtering out the significant markers from GWAS
stp_gwas_markers <- unlist(as.vector(stp_gwas[,1]))
sol_Mo_gwas_markers <- unlist(as.vector(sol_Mo_gwas[,1]))
PNZ_gwas_markers <- unlist(as.vector(PNZ_gwas[,1]))
POL_gwas_markers <- unlist(as.vector(POL_gwas[,1]))

##### ASSUMPTION - the missing makers are homozygous
# Subsetting markers - replacing -9 with 0
stp_markers <- SNP_markers[stp_gwas_markers]
stp_markers <- replace(stp_markers, stp_markers == -9,0)

sol_Mo_markers <- SNP_markers[sol_Mo_gwas_markers]
sol_Mo_markers <- replace(sol_Mo_markers, sol_Mo_markers == -9,0)

PNZ_markers <- SNP_markers[PNZ_gwas_markers]
PNZ_markers <- replace(PNZ_markers, PNZ_markers == -9,0)

POL_markers <- SNP_markers[POL_gwas_markers]
POL_markers <- replace(POL_markers, POL_markers == -9,0)


##### Subsetting Testing and Training data for markers
set.seed(123)

tot_train <- as.matrix(sample(1:nrow(tot_markers), 1035))
tot_test <- setdiff(1:nrow(tot_markers), tot_train)

lab_train <- as.matrix(sample(1:nrow(lab_markers), 1035))
lab_test <- setdiff(1:nrow(lab_markers), lab_train)

sol_Mo_train <- as.matrix(sample(1:nrow(sol_Mo_markers), 1035))
sol_Mo_test <- setdiff(1:nrow(sol_Mo_markers), sol_Mo_train)


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


sol_Mo_train_phenotype <- as.matrix(sol_Mo[sol_Mo_train, ])
sol_Mo_train_marker <- as.matrix(sol_Mo_markers[sol_Mo_train, ], K = NULL)

# Removing NAs from training and putting them into testing
#nas <- which(is.na(sol_Mo_train_phenotype))

#sol_Mo_train_phenotype <- as.matrix(sol_Mo_train_phenotype[-nas, ])
#sol_Mo_train_marker <- as.matrix(sol_Mo_train_marker[-nas, ])
#sol_Mo_train <- sol_Mo_train[-nas,]
#sol_Mo_test <- sort(c(sol_Mo_test, nas))
sol_Mo_test_phenotype <- as.matrix(sol_Mo[sol_Mo_test, ])
sol_Mo_test_marker <- as.matrix(sol_Mo_markers[sol_Mo_test, ], K = NULL)





# Training the model
model_tot <- mixed.solve(y = tot_train_phenotype, Z = tot_train_marker)

model_lab <- mixed.solve(y = lab_train_phenotype, Z = lab_train_marker)

model_sol_Mo <- mixed.solve(y = sol_Mo_train_phenotype, Z = sol_Mo_train_marker)


# Marker effects
marker_effects_tot <- as.matrix(model_tot$u)

marker_effects_lab <- as.matrix(model_lab$u)

marker_effects_sol_Mo <- as.matrix(model_sol_Mo$u)


# Running BLUE
BLUE_tot <- as.vector(model_tot$beta)

BLUE_lab <- as.vector(model_lab$beta)

BLUE_sol_Mo <- as.vector(model_sol_Mo$beta)


# Prediction
predict_train_tot <- as.matrix(tot_train_marker) %*% marker_effects_tot
predict_test_tot <- as.matrix(tot_test_marker) %*% marker_effects_tot
predict_train_result_tot <- as.vector((predict_train_tot[,1])+BLUE_tot)
predict_test_result_tot <- as.vector((predict_test_tot[,1])+BLUE_tot)

predict_train_lab <- as.matrix(lab_train_marker) %*% marker_effects_lab
predict_test_lab <- as.matrix(lab_test_marker) %*% marker_effects_lab
predict_train_result_lab <- as.vector((predict_train_lab[,1])+BLUE_lab)
predict_test_result_lab <- as.vector((predict_test_lab[,1])+BLUE_lab)

predict_train_sol_Mo <- as.matrix(sol_Mo_train_marker) %*% marker_effects_sol_Mo
predict_test_sol_Mo <- as.matrix(sol_Mo_test_marker) %*% marker_effects_sol_Mo
predict_train_result_sol_Mo <- as.vector((predict_train_sol_Mo[,1])+BLUE_sol_Mo)
predict_test_result_sol_Mo <- as.vector((predict_test_sol_Mo[,1])+BLUE_sol_Mo)


# Correlation
corr_tot_train <- cor.test(as.vector(tot_train_phenotype), predict_train_result_tot, use = "spearman")
corr_tot_test <- cor(as.vector(tot_test_phenotype), predict_test_result_tot, use = "complete")

corr_lab_train <- cor.test(as.vector(lab_train_phenotype), predict_train_result_lab, use = "spearman")
corr_lab_test <- cor(as.vector(lab_test_phenotype), predict_test_result_lab, use = "complete")

corr_sol_Mo_train <- cor.test(as.vector(sol_Mo_train_phenotype), predict_train_result_sol_Mo, use = "spearman")
corr_sol_Mo_test <- cor(as.vector(sol_Mo_test_phenotype), predict_test_result_sol_Mo, use = "complete")


# Visualize
plot(tot_test_phenotype, predict_test_result_tot)
abline(lm(tot_test_phenotype~predict_test_result_tot), col = "red")

plot(lab_test_phenotype, predict_test_result_lab)
plot(sol_Mo_test_phenotype, predict_test_result_sol_Mo)
abline(lm(sol_Mo_test_phenotype~ predict_test_result_sol_Mo), col = "red")

# Visualizing side by side
tot_train_test <- cbind(tot_test_phenotype, predict_test_result_tot)
lab_train_test <- cbind(lab_test_phenotype, predict_test_result_lab)
sol_Mo_train_test <- cbind(sol_Mo_test_phenotype, predict_test_result_sol_Mo)

x <- cbind(sol_Mo_train_test,sol_Mo_train_test_DNN)
x <- x[,-3]
colnames(x) <- c("soilP","rrBLUP","CNN")
x <- as.data.frame(x)
str(x)
x$CNN <- x$CNN * 2





# SNP interaction data 
snp_interactions is the matrix of interactions
A <- A.mat(snp_interactions, model = "Additive")
result <- mixed.solve(y = y, X = X, K = A)

