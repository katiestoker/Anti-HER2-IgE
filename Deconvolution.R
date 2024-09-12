
## Author: Roman Laddach
# A script to deconvolute bulk TPM data and visualize the levels using a heatmap.  

# Set up environment 
library("dplyr") 
library("reshape2") 
library("ggplot2") 
library("pheatmap") 
library("ConsensusTME") 

# Load rat information
rat_info = read.csv("/path_to/rat_information.csv")

# Load TPM counts 
TPM_counts = readRDS("/path_to/TPM_counts_collapsed_multiplicates_human_genes.RDS")

# Adding rownames to TPM_counts and aligning the rat name to match the rat_info  
rownames(TPM_counts) = TPM_counts$human_symbol
TPM_counts$human_symbol = NULL
TPM_counts$symbol = NULL

# Preprocessing
rat_info$Animal_ID
rat_info$Animal_ID = gsub("-", "", rat_info$Animal_ID)
rownames(rat_info) = rat_info$Animal_ID
rat_info$Animal_ID

colnames(TPM_counts)
colnames(TPM_counts) = gsub("\\.", "", colnames(TPM_counts))

# Check if the order is identical  
identical(rownames(rat_info), colnames(TPM_counts))

# Rearrange 
TPM_counts = TPM_counts[, rat_info$Animal_ID]
identical(rownames(rat_info), colnames(TPM_counts))

## Deconvolution
bulkExpMatrix = as.matrix(TPM_counts)
bulkExpMatrix[1:5,1:5]
all_samples_deconvolution = consensusTMEAnalysis(bulkExpMatrix, cancer = "BRCA", statMethod = "ssgsea")

# Results 
all_samples_deconvolution[1:5,1:5]

# Transpose the results 
all_samples_deconvolution_t = as.data.frame(t(all_samples_deconvolution))
head(all_samples_deconvolution_t)