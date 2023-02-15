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


###############
#Importing Data
###############
data_tbl <- read.csv("only_complete_Clinical_Data_Patient_Level-TCells.csv", header = T)

# Just the complete data
head(data_tbl)
names(data_tbl)
table(data_tbl$BestDiseaseResponse)
table(data_tbl$Transfection_Method)


data_tbl <- data_tbl %>%
  dplyr::mutate(binary_response = case_when(BestDiseaseResponse == "CR" ~ 1, TRUE ~ 0))

data_tbl <- data_tbl[,-7]

# Finding unique variables and converting them to numbers. 
for (i in 1:ncol(data_tbl)) {
  if (is.character(data_tbl[,i])) {
    data_tbl[,i] <- factor(data_tbl[,i])
    data_tbl[,i] <- as.numeric(data_tbl[,i])
  }
}

data_tbl$Patient_Age <- (data_tbl$Patient_Age-min(data_tbl$Patient_Age))/(max(data_tbl$Patient_Age)-min(data_tbl$Patient_Age))
ggplot(data_tbl, aes(x=data_tbl$Patient_Malignancy,y=data_tbl$Patient_Age)) + geom_point()

cor.test(data_tbl$Patient_Malignancy,data_tbl$Patient_Age)



ggplot(data_tbl, aes(x=data_tbl$Patient_Malignancy,y=data_tbl$Patient_Age)) + geom_point()

head(data_tbl)
str(data_tbl)

#Modeling
fit_xgboost <- boost_tree(learn_rate = 0.3) %>%
  set_mode("regression") %>%
  set_engine("xgboost") %>%
  fit(binary_response ~. , data = data_tbl)
fit_xgboost

hist(data_tbl$Patient_Age)

plot()
# Explainer
explainer <- DALEX::explain(
  model = fit_xgboost,
  data = data_tbl,
  y = data_tbl$binary_response,
  label = "XGBoost"
)


# MODEL STUDIO visualization
modelStudio::modelStudio(explainer)

table(data_tbl$Patient_Malignancy)
unique(data_tbl$Patient_Malignancy)
