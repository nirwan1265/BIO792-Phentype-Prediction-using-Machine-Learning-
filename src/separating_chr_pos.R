# Change working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Github/BIO792-Phentype-Prediction-using-Machine-Learning-/data/nirwan_data")
# Phenotypes
phenotypes <- read.csv("Sorghum_allphospho_africa.csv")
# Total phosphorus
tot <- phenotypes["tot"]
#hist(tot[,1])

# Solubility
lab <- phenotypes["lab"]
#hist(lab[,1])

# Phosphorus retention
sol_VL <- phenotypes["sol_VL"]
#hist(sol_VL[,1])


# Loading GWAS files and only selecting the significant snps
tot_gwas <- vroom("tot_LMM.txt") %>% select(rs,p_wald) %>% filter(p_wald <= 0.05) %>% select(1)
lab_gwas <- vroom("lab_LMM.txt") %>% select(rs,p_wald) %>% filter(p_wald <= 0.05) %>% select(1)
sol_VL_gwas <- vroom("sol_VL_LMM.txt") %>% select(rs,p_wald) %>% filter(p_wald <= 0.05) %>% select(1)

# extract chromosome number and position
sol_VL_gwas$chr <- gsub(pattern = "^S(\\d+)_.*", replacement = "Chr0\\1", x = sol_VL_gwas$rs)
sol_VL_gwas$chr <- gsub(pattern = "^S(\\d{2})_.*", replacement = "Chr\\1", x = sol_VL_gwas$chr)
sol_VL_gwas$chr <- gsub(pattern = "^Chr0(\\d{2})", replacement = "Chr\\1", x = sol_VL_gwas$chr)

sol_VL_gwas$pos <- gsub(pattern = "^S\\d+_(\\d+)$", replacement = "\\1", x = sol_VL_gwas$rs)

sol_VL_gwas <- sol_VL_gwas[,-1]


# get unique chromosome names
chr_names <- unique(sol_VL_gwas$chr)

# loop over chromosome names and assign to separate variables
for (chr_name in chr_names) {
  # subset data frame for the current chromosome
  chr_df <- sol_VL_gwas[sol_VL_gwas$chr == chr_name, c("chr", "pos")]
  
  # construct file name
  file_name <- paste0("sol_VL_", chr_name, "_gwas.txt")
  
  # write to file
  write.table(chr_df, file = file_name, row.names = FALSE, quote = FALSE)
}
