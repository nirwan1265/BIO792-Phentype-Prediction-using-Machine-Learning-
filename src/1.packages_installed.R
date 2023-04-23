library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(caret)
library(randomForest)
library(rrBLUP)
library(tensorflow)
library(keras)
library(kerasR)
library(doParallel)
library(foreach)
library(vroom)
library(dplyr)
library(reticulate)
library(e1071)
library(glmnet)
library(gbm)
library(xgboost)

#install.packages("BayesR", repos="http://R-Forge.R-project.org")
#path_to_python <- install_python()
#Note that if you already have Python installed, you don’t need to call install_python() and instead can just supply an absolute path to the Python executable.
virtualenv_create("r-reticulate", python = path_to_python)

#install.packages("tensorflow")


#library(tensorflow)
#install_tensorflow(envname = "r-reticulate")
#install.packages("keras")
#library(keras)
#install_keras()
#install.packages("kerasR")
#devtools::install_github("rstudio/keras")


