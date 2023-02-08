library(glmnet)
library(modelStudio)
library(DALEX)
library(tidyverse)
library(tidymodels)
library(xgboost)
library(dplyr)


##############
#XGBOOST
##############
setwd(paste0(getwd(),"/data/lucas_data"))
system("ls")

###############
#Importing Data
###############
data_tbl <- read.csv("only_complete_Clinical_Data_Patient_Level-TCells.csv", header = T)
# Just the complete data
head(data_tbl)
names(data_tbl)


# Finding unique variables and converting them to numbers. 
for (i in 1:ncol(data_tbl)) {
  if (is.character(data_tbl[,i])) {
    data_tbl[,i] <- factor(data_tbl[,i])
    data_tbl[,i] <- as.numeric(data_tbl[,i])
  }
}
head(data_tbl)
names(data_tbl)

#Modeling
fit_xgboost <- boost_tree(learn_rate = 0.3, trees = 10) %>%
  set_mode("regression") %>%
  set_engine("xgboost") %>%
  fit(BestDiseaseResponse ~. , data = data_tbl)
fit_xgboost


# Explainer
explainer <- DALEX::explain(
  model = fit_xgboost,
  data = data_tbl,
  y = data_tbl$BestDiseaseResponse,
  label = "XGBoost"
)


# MODEL STUDIO visualization
modelStudio::modelStudio(explainer)

