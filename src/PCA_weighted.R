# Library packages
library(pls)
library(caret)
library(pcaL1)
library(vroom)
library(glmnet)
library(modelStudio)
library(DALEX)
library(tidyverse)
library(tidymodels)
library(xgboost)

# Setting working directory for data
setwd("./data/nirwan_data")

# Phenotype data
pheno <- read.csv("Sorghum_allphospho_africa.csv")

# Loading PCA file
GPCA <- readRDS("pca_sorghum.RDS")

#Make a table of eigen values
j <- 1
for(i in paste0("GPCA")){
  d = get(i)
  assign(paste0("EV",sprintf("%02d",j)), data.frame(sample.id = d$sample.id,
                                                     EV1 = d$eigenvect[,1],
                                                     EV2 = d$eigenvect[,2],
                                                     EV3 = d$eigenvect[,3],
                                                     EV4 = d$eigenvect[,4],
                                                     EV5 = d$eigenvect[,5],
                                                     EV6 = d$eigenvect[,6],
                                                     EV7 = d$eigenvect[,7],
                                                     EV8 = d$eigenvect[,8],
                                                     EV9 = d$eigenvect[,9],
                                                     EV10 = d$eigenvect[,10],
                                                     EV11 = d$eigenvect[,11],
                                                     EV12 = d$eigenvect[,12],
                                                     EV13 = d$eigenvect[,13],
                                                     EV14 = d$eigenvect[,14],
                                                     EV15 = d$eigenvect[,15],
                                                     EV16 = d$eigenvect[,16],
                                                     EV17 = d$eigenvect[,17],
                                                     EV18 = d$eigenvect[,18],
                                                     EV19 = d$eigenvect[,19],
                                                     EV20 = d$eigenvect[,20],
                                                     EV21 = d$eigenvect[,21],
                                                     EV22 = d$eigenvect[,22],
                                                     EV23 = d$eigenvect[,23],
                                                     EV24 = d$eigenvect[,24],
                                                     EV25 = d$eigenvect[,25],
                                                     stringsAsFactors = FALSE))
  assign(i,d)
  j = j + 1
}


# Variance Proportion (%)
pc.percent <- GPCA$varprop*100
pc.percent[1:25]
round(pc.percent)[1:25]


# Top 13 EV explains atleast some variance of the Genotype
EV01 <- EV01[,-c(1,14:26)]
EV01 <- EV01[,2]


# Giving weights to the Eigen Vector:
W = matrix(c(2,2,1,1,1,1,1,1,1,1,1,1), ncol = 1)


#Combining Eigen Vector with weight 
composite_pca <- as.matrix(EV01)
composite_pca <- as.data.frame(rowSums(as.matrix(composite_pca) %*% W))
head(composite_pca)


##############
# Partial Least Square
##############

# Fit a PLS model using the PCs as the independent variables and the phenotypes as the dependent variables
pls_model <- plsr(as.matrix(composite_pca)~as.matrix(pheno[,2:29]), ncomp = 1)
summary(pls_model)

# Print the variable importance in the projection (VIP) scores for the phenotypes
imp <- as.data.frame(varImp(pls_model))
imp_df <- as.data.frame(imp)
imp_df$pheno <- rownames(imp)
imp_sorted <- imp_df[order(-imp_df$Overall),]
rownames(imp_sorted) <- NULL
imp_sorted$pheno <- sapply(strsplit(imp_sorted$pheno, ")"), "[",2)
imp_sorted






##############
#XGBOOST
##############
data_tbl <- as_tibble(cbind(as.data.frame(composite_pca[,1]), pheno[,2:29]))
names(data_tbl)[1] <- "pca"


#Model
fit_xgboost <- boost_tree(learn_rate = 0.3) %>%
  set_mode("regression") %>%
  set_engine("xgboost") %>%
  fit(pca ~. , data = data_tbl)

fit_xgboost


# Explainer
explainer <- DALEX::explain(
  model = fit_xgboost,
  data = data_tbl,
  y = data_tbl$pca,
  label = "XGBoost"
)


# MODEL STUDIO

modelStudio::modelStudio(explainer)
