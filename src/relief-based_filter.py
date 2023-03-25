# Required Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import mutual_info_regression, mutual_info_classif
from sklearn.preprocessing import LabelEncoder
import os

lab_markers_path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/selected_SNPs/lab_markers.txt"
sol_VL_markers_path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/selected_SNPs/sol_VL_markers.txt"
tot_markers_path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/selected_SNPs/tot_markers.txt"
pheno_path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/Sorghum_allphospho_africa.csv"

lab_markers = pd.read_csv(lab_markers_path, sep = " ")
sol_VL_markers = pd.read_csv(sol_VL_markers_path, sep = " ")
tot_markers = pd.read_csv(tot_markers_path, sep = " ")
pheno = pd.read_csv(pheno_path)

lab_pheno = pheno["lab"]
sol_VL_pheno = pheno["sol_VL"]
tot_pheno = pheno["tot"]


type(lab_markers)

lab = pd.concat([lab_pheno,lab_markers], axis=1)
sol_VL = pd.concat([sol_VL_pheno,sol_VL_markers], axis=1)
tot = pd.concat([tot_pheno,tot_markers], axis=1)


# Columns 
lab_col = lab.shape[1]
sol_VL_col = sol_VL.shape[1]
tot_col = tot.shape[1]


#Creating the input features X and target variable y
X_lab = lab.iloc[:,1:lab_col]
y_lab = lab.iloc[:,0].values

X_sol_VL = sol_VL.iloc[:,1:sol_VL_col]
y_sol_VL = sol_VL.iloc[:,0].values

X_tot = tot.iloc[:,1:tot_col]
y_tot = tot.iloc[:,0].values

# Create a data set with all the input features after converting them to numeric including target variable
full_data_tot= X_tot.copy()
full_data_tot['tot']= y_tot
full_data_tot.head(8)

full_data_sol_VL= X_sol_VL.copy()
full_data_sol_VL['sol_VL']= y_sol_VL
full_data_sol_VL.head(8)

full_data_lab= X_lab.copy()
full_data_lab['lab']= y_lab
full_data_lab.head(8)


# Applying step 1 of the filter method
# Identify input features having high correlation with target variable

importances_lab = full_data_lab.drop("lab", axis=1).apply(lambda x: x.corr(full_data_lab.lab))
indices_lab = np.argsort(importances_lab)
print(importances_lab[indices_lab])

importances_tot = full_data_tot.drop("tot", axis=1).apply(lambda x: x.corr(full_data_tot.tot))
indices_tot = np.argsort(importances_tot)
print(importances_tot[indices_tot])

importances_sol_VL = full_data_sol_VL.drop("sol_VL", axis=1).apply(lambda x: x.corr(full_data_sol_VL.sol_VL))
indices_sol_VL = np.argsort(importances_sol_VL)
print(importances_sol_VL[indices_sol_VL])


#Plotting 
exclude_lab = ['lab']
names_lab = [col for col in lab.columns if col not in exclude_lab]

exclude_tot = ['tot']
names_tot = [col for col in tot.columns if col not in exclude_tot]

exclude_sol_VL = ['sol_VL']
names_sol_VL = [col for col in sol_VL.columns if col not in exclude_sol_VL]


# We want to keep features with only a high correlation with the target variable. This implies that the input feature has a high influence in predicting the target variable.
# We set the threshold to the absolute value of 0.3. We keep input features only if the correlation of the input feature with the target variable is greater than 0.4

X_lab = lab[[names_lab[i] for i in range(len(indices_lab)) if np.abs(importances_lab[i])>=0.2]]

X_tot = tot[[names_tot[i] for i in range(len(indices_tot)) if np.abs(importances_tot[i])>=0.2]]

X_sol_VL = sol_VL[[names_sol_VL[i] for i in range(len(indices_sol_VL)) if np.abs(importances_sol_VL[i])>=0.2]]


# Saving
X_lab.to_csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/selected_SNPs/filtered_lab_markers.txt", sep="\t", index=False)
X_tot.to_csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/selected_SNPs/filtered_tot_markers.txt", sep="\t", index=False)
X_sol_VL.to_csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data/selected_SNPs/filtered_sol_VL_markers.txt", sep="\t", index=False)

