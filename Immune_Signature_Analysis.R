# A script to investigate immune cell enrichment in treated/non-treated rats.

# Set-up environment 
library(ConsensusTME)
library(dplyr)
library(tidyverse)
sessionInfo()

# View the consensus gene sets available for different cancers 
data(consensusGeneSets)

# Load the breast cancer consensus gene sets for the immune cells of interest

# M1
m1_cells = consensusGeneSets[["BRCA"]][["Macrophages_M1"]]

# M2
m2_cells = consensusGeneSets[["BRCA"]][["Macrophages_M2"]]

# CD8 T Cells
t_cells = consensusGeneSets[["BRCA"]][["T_cells_CD8"]]

# CD4 T cell
CD4_cells = consensusGeneSets[["BRCA"]][["T_cells_CD4"]]

# NK_Cells
NK_cells = consensusGeneSets[["BRCA"]][["NK_cells"]]

# Cytotoxic T cell
cytotoxic_cells = consensusGeneSets[["BRCA"]][["Cytotoxic_cells"]]

# Gamma/Delta
gamma_delta_cells = consensusGeneSets[["BRCA"]][["T_cells_gamma_delta"]]


# Define the pheatmap.scale function

pheatmap.scale <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  return((x - m) / s)
}


# Load the required datasets 

# Load the DEGS 
DEGs_AbC_vs_PBS = read.csv("/path_to/DEGs_AbC_vs_PBS.csv", row.names = 1)
head(DEGs_AbC_vs_PBS)

# Load the metadata 
rat_info = read.csv("/path_to/rat_information.csv")

# Subset rat info for only pbs vs AbC data 
rat_info_sub = rat_info %>% filter(Treatment_short %in% c("PBS", "Ab_C"))

# Load TPM count
TPM = readRDS("/path_to/TPM_counts_collapsed_multiplicates_human_genes.RDS")
rownames(TPM) = TPM$human_symbol
TPM$human_symbol = NULL
TPM$symbol = NULL

# Change the column names in TPM data to match animal ID in rat info metadata (dash not full stop)
colnames(TPM) = gsub("\\.", "-", colnames(TPM))
TPM_sub = TPM[, rat_info_sub$Animal_ID]


# Run code to interrogate each cell type individually 

## M1 ##
m1_cells_cell_expression = rownames(DEGs_AbC_vs_PBS)[rownames(DEGs_AbC_vs_PBS) %in% m1_cells]

# Subset the TPM data to include only GOI (monocyte signature expression) 
TPM_sub_m1_cells = TPM_sub[m1_cells_cell_expression, ]

# Determine heatmap annotation
rat_info_sub %>% group_by(Treatment_short)
sample_annotation = rat_info_sub[("Treatment_short")]
rownames(sample_annotation) <- colnames(TPM_sub_m1_cells)
TPM_sub_m1_cells = TPM_sub_m1_cells[do.call(order, c(TPM_sub_m1_cells, list(decreasing=TRUE))),]

# Order the TPM/DEG selected data from high to low
TPM_sub_m1 = TPM_sub_m1_cells[do.call(order, c(TPM_sub_m1_cells, list(decreasing=TRUE))),]

# Merge the DEGS/TPM data with DEGs metadata
TPM_sub_m1_merged <- merge(TPM_sub_m1, DEGs_AbC_vs_PBS, by = "row.names")

# Check which degs are significant
significant_degs = TPM_sub_m1_merged %>% filter(padj<=0.05) %>% arrange(desc(log2FoldChange))
significant_degs = column_to_rownames(significant_degs, var = "Row.names")

# Append dataframne to show DEG is significant
TPM_sub_m1_merged$significance = " "
TPM_sub_m1_merged = column_to_rownames(TPM_sub_m1_merged, var = "Row.names")
TPM_sub_m1_merged$significance <- ifelse(rownames(TPM_sub_m1_merged) %in% rownames(significant_degs), "*", " ")

# Change rownames to first column
TPM_sub_m1_merged <- rownames_to_column(TPM_sub_m1_merged, "genes") # Apply rownames_to_column

# Make a new column combining gene name and significance 
TPM_sub_m1_merged$row_names_combined <- paste0(TPM_sub_m1_merged$genes, " ", TPM_sub_m1_merged$significance)

# Set the new column as row names and remove unnecessary columns
rownames(TPM_sub_m1_merged) <- TPM_sub_m1_merged$row_names_combined
TPM_sub_m1_merged <- TPM_sub_m1_merged[, -c(1, ncol(TPM_sub_m1_merged))]

TPM_sub_m1_merged_heatmap <- TPM_sub_m1_merged[, 1:13]

# Manually Scale the Data 
TPM_sub_M1_merged_heatmap.scaled <- pheatmap.scale(TPM_sub_m1_merged_heatmap)

# Export the data 
write.csv(TPM_sub_M1_merged_heatmap.scaled, "/path_to/TPM_sub_M1_merged_heatmap.scaled.csv")


## M2 ##
m2_cells_cell_expression = rownames(DEGs_AbC_vs_PBS)[rownames(DEGs_AbC_vs_PBS) %in% m2_cells]

# Subset the TPM data to include only GOI (monocyte signature expression) 
TPM_sub_m2_cells = TPM_sub[m2_cells_cell_expression, ]
TPM_sub_m2_cells = TPM_sub_m2_cells %>% filter(rowSums(across(where(is.numeric)))!=0)

# Determine heatmap annotation
rat_info_sub %>% group_by(Treatment_short)
sample_annotation = rat_info_sub[("Treatment_short")]
rownames(sample_annotation) <- colnames(TPM_sub_m2_cells)
TPM_sub_m2_cells = TPM_sub_m2_cells[do.call(order, c(TPM_sub_m2_cells, list(decreasing=TRUE))),]

# Order the TPM/DEG selected data from high to low
TPM_sub_m2 = TPM_sub_m2_cells[do.call(order, c(TPM_sub_m2_cells, list(decreasing=TRUE))),]

# Merge the DEGS/TPM data with DEGs metadata
TPM_sub_m2_merged <- merge(TPM_sub_m2, DEGs_AbC_vs_PBS, by = "row.names")

# Check which degs are significant - 13 DEGS
significant_degs = TPM_sub_m2_merged %>% filter(padj<=0.05) %>% arrange(desc(log2FoldChange))
significant_degs = column_to_rownames(significant_degs, var = "Row.names")

# Append dataframne to show DEG is significant
TPM_sub_m2_merged$significance = " "
TPM_sub_m2_merged = column_to_rownames(TPM_sub_m2_merged, var = "Row.names")
TPM_sub_m2_merged$significance <- ifelse(rownames(TPM_sub_m2_merged) %in% rownames(significant_degs), "*", " ")

# Change rownames to first column
TPM_sub_m2_merged <- rownames_to_column(TPM_sub_m2_merged, "genes") # Apply rownames_to_column

# Make a new column combining gene name and significance 
TPM_sub_m2_merged$row_names_combined <- paste0(TPM_sub_m2_merged$genes, " ", TPM_sub_m2_merged$significance)

# Set the new column as row names and remove unnecessary columns
rownames(TPM_sub_m2_merged) <- TPM_sub_m2_merged$row_names_combined
TPM_sub_m2_merged <- TPM_sub_m2_merged[, -c(1, ncol(TPM_sub_m2_merged))]
TPM_sub_m2_merged_heatmap <- TPM_sub_m2_merged[, 1:13]

# Manually Scale the Data 
TPM_sub_M2_merged_heatmap.scaled <- pheatmap.scale(TPM_sub_m2_merged_heatmap)

# Export the data 
write.csv(TPM_sub_M2_merged_heatmap.scaled, "/path_to/TPM_sub_M2_merged_heatmap_scaled.csv")


## CD8 T Cells ##
T_cells_cell_expression = rownames(DEGs_AbC_vs_PBS)[rownames(DEGs_AbC_vs_PBS) %in% t_cells]

# Subset the TPM data to include only GOI 
TPM_sub_T_cells = TPM_sub[T_cells_cell_expression, ]
TPM_sub_T_cells = TPM_sub_T_cells %>% filter(rowSums(across(where(is.numeric)))!=0)

# Determine heatmap annotation
rat_info_sub %>% group_by(Treatment_short)
sample_annotation = rat_info_sub[("Treatment_short")]
rownames(sample_annotation) <- colnames(TPM_sub_T_cells)
TPM_sub_T_cells = TPM_sub_T_cells[do.call(order, c(TPM_sub_T_cells, list(decreasing=TRUE))),]

# Order the TPM/DEG selected data from high to low
TPM_sub_T = TPM_sub_T_cells[do.call(order, c(TPM_sub_T_cells, list(decreasing=TRUE))),]

# Merge the DEGS/TPM data with DEGs metadata
TPM_sub_T_merged <- merge(TPM_sub_T, DEGs_AbC_vs_PBS, by = "row.names")

# Check which degs are significant - 13 DEGS
significant_degs = TPM_sub_T_merged %>% filter(padj<=0.05) %>% arrange(desc(log2FoldChange))
significant_degs = column_to_rownames(significant_degs, var = "Row.names")

# Append dataframne to show DEG is significant
TPM_sub_T_merged$significance = " "
TPM_sub_T_merged = column_to_rownames(TPM_sub_T_merged, var = "Row.names")
TPM_sub_T_merged$significance <- ifelse(rownames(TPM_sub_T_merged) %in% rownames(significant_degs), "*", " ")

# Change rownames to first column
TPM_sub_T_merged <- rownames_to_column(TPM_sub_T_merged, "genes") # Apply rownames_to_column

# Make a new column combining gene name and significance 
TPM_sub_T_merged$row_names_combined <- paste0(TPM_sub_T_merged$genes, " ", TPM_sub_T_merged$significance)

# Set the new column as row names and remove unnecessary columns
rownames(TPM_sub_T_merged) <- TPM_sub_T_merged$row_names_combined
TPM_sub_T_merged <- TPM_sub_T_merged[, -c(1, ncol(TPM_sub_T_merged))]
TPM_sub_T_merged_heatmap <- TPM_sub_T_merged[, 1:13]

# Manually Scale the Data 
TPM_sub_T_merged_heatmap.scaled <- pheatmap.scale(TPM_sub_T_merged_heatmap)

# Export the data 
write.csv(TPM_sub_T_merged_heatmap.scaled, "/path_to/TPM_sub_T_cell_merged_heatmap_scaled.csv")


## NK Cell ##
NK_cells_cell_expression = rownames(DEGs_AbC_vs_PBS)[rownames(DEGs_AbC_vs_PBS) %in% NK_cells]

# Subset the TPM data to include only GOI 
TPM_sub_NK_cells = TPM_sub[NK_cells_cell_expression, ]
TPM_sub_NK_cells = TPM_sub_NK_cells %>% filter(rowSums(across(where(is.numeric)))!=0)

# Determine heatmap annotation
rat_info_sub %>% group_by(Treatment_short)
sample_annotation = rat_info_sub[("Treatment_short")]
rownames(sample_annotation) <- colnames(TPM_sub_NK_cells)
TPM_sub_NK_cells = TPM_sub_NK_cells[do.call(order, c(TPM_sub_NK_cells, list(decreasing=TRUE))),]

# Order the TPM/DEG selected data from high to low
TPM_sub_NK_cells = TPM_sub_NK_cells[do.call(order, c(TPM_sub_NK_cells, list(decreasing=TRUE))),]

# Merge the DEGS/TPM data with DEGs metadata
TPM_sub_NK_cells_merged <- merge(TPM_sub_NK_cells, DEGs_AbC_vs_PBS, by = "row.names")

# Check which degs are significant 
significant_degs = TPM_sub_NK_cells_merged %>% filter(padj<=0.05) %>% arrange(desc(log2FoldChange))
significant_degs = column_to_rownames(significant_degs, var = "Row.names")

# Append dataframne to show DEG is significant
TPM_sub_B_cell_merged$significance = " "
TPM_sub_NK_cells_merged = column_to_rownames(TPM_sub_NK_cells_merged, var = "Row.names")
TPM_sub_NK_cells_merged$significance <- ifelse(rownames(TPM_sub_NK_cells_merged) %in% rownames(significant_degs), "*", " ")

# Change rownames to first column
TPM_sub_NK_cells_merged <- rownames_to_column(TPM_sub_NK_cells_merged, "genes") # Apply rownames_to_column

# Make a new column combining gene name and significance 
TPM_sub_NK_cells_merged$row_names_combined <- paste0(TPM_sub_NK_cells_merged$genes, " ", TPM_sub_NK_cells_merged$significance)

# Set the new column as row names and remove unnecessary columns
rownames(TPM_sub_NK_cells_merged) <- TPM_sub_NK_cells_merged$row_names_combined
TPM_sub_NK_cells_merged <- TPM_sub_NK_cells_merged[, -c(1, ncol(TPM_sub_NK_cells_merged))]
TPM_sub_NK_cells_merged_heatmap <- TPM_sub_NK_cells_merged[, 1:13]

# Manually Scale the Data 
TPM_sub_NK_cell_merged_heatmap.scaled <- pheatmap.scale(TPM_sub_NK_cells_merged_heatmap)

# Export the data 
write.csv(TPM_sub_NK_cell_merged_heatmap.scaled, "/path_to/TPM_sub_NK_cell_merged_heatmap_scaled.csv")



##  CD4 Cell ##
CD4_cells_cell_expression = rownames(DEGs_AbC_vs_PBS)[rownames(DEGs_AbC_vs_PBS) %in% CD4_cells]

# Subset the TPM data to include only GOI 
TPM_sub_CD4_cells = TPM_sub[CD4_cells_cell_expression, ]
TPM_sub_CD4_cells = TPM_sub_CD4_cells %>% filter(rowSums(across(where(is.numeric)))!=0)

# Determine heatmap annotation
rat_info_sub %>% group_by(Treatment_short)
sample_annotation = rat_info_sub[("Treatment_short")]
rownames(sample_annotation) <- colnames(TPM_sub_CD4_cells)
TPM_sub_CD4_cells = TPM_sub_CD4_cells[do.call(order, c(TPM_sub_CD4_cells, list(decreasing=TRUE))),]

# Order the TPM/DEG selected data from high to low
TPM_sub_CD4_cells = TPM_sub_CD4_cells[do.call(order, c(TPM_sub_CD4_cells, list(decreasing=TRUE))),]

# Merge the DEGS/TPM data with DEGs metadata
TPM_sub_CD4_cells_merged <- merge(TPM_sub_CD4_cells, DEGs_AbC_vs_PBS, by = "row.names")

# Check which degs are significant 
significant_degs = TPM_sub_CD4_cells_merged %>% filter(padj<=0.05) %>% arrange(desc(log2FoldChange))
significant_degs = column_to_rownames(significant_degs, var = "Row.names")

# Append dataframne to show DEG is significant
TPM_sub_CD4_cells_merged$significance = " "
TPM_sub_CD4_cells_merged = column_to_rownames(TPM_sub_CD4_cells_merged, var = "Row.names")
TPM_sub_CD4_cells_merged$significance <- ifelse(rownames(TPM_sub_CD4_cells_merged) %in% rownames(significant_degs), "*", " ")

# Change rownames to first column
TPM_sub_CD4_cells_merged <- rownames_to_column(TPM_sub_CD4_cells_merged, "genes") # Apply rownames_to_column

# Make a new column combining gene name and significance 
TPM_sub_CD4_cells_merged$row_names_combined <- paste0(TPM_sub_CD4_cells_merged$genes, " ", TPM_sub_CD4_cells_merged$significance)

# Set the new column as row names and remove unnecessary columns
rownames(TPM_sub_CD4_cells_merged) <- TPM_sub_CD4_cells_merged$row_names_combined
TPM_sub_CD4_cells_merged <- TPM_sub_CD4_cells_merged[, -c(1, ncol(TPM_sub_CD4_cells_merged))]
TPM_sub_CD4_cells_merged_heatmap <- TPM_sub_CD4_cells_merged[, 1:13]

# Manually Scale the Data 
TPM_sub_CD4_cells_merged_heatmap.scaled <- pheatmap.scale(TPM_sub_CD4_cells_merged_heatmap)

# Export the data 
write.csv(TPM_sub_CD4_cells_merged_heatmap.scaled, "/path_to/TPM_sub_CD4_cell_merged_heatmap_scaled.csv")


##  Cytotoxic T Cell ##
cyto_cells_cell_expression = rownames(DEGs_AbC_vs_PBS)[rownames(DEGs_AbC_vs_PBS) %in% cytotoxic_cells]

# Subset the TPM data to include only GOI 
TPM_sub_cyto_cells = TPM_sub[cyto_cells_cell_expression, ]
TPM_sub_cyto_cells = TPM_sub_cyto_cells %>% filter(rowSums(across(where(is.numeric)))!=0)

# Determine heatmap annotation
rat_info_sub %>% group_by(Treatment_short)
sample_annotation = rat_info_sub[("Treatment_short")]
rownames(sample_annotation) <- colnames(TPM_sub_cyto_cells)
TPM_sub_cyto_cells = TPM_sub_cyto_cells[do.call(order, c(TPM_sub_cyto_cells, list(decreasing=TRUE))),]

# Order the TPM/DEG selected data from high to low
TPM_sub_cyto_cells = TPM_sub_cyto_cells[do.call(order, c(TPM_sub_cyto_cells, list(decreasing=TRUE))),]
TPM_sub_cyto_cells[1:10, ]

# Merge the DEGS/TPM data with DEGs metadata
TPM_sub_cyto_cells_merged <- merge(TPM_sub_cyto_cells, DEGs_AbC_vs_PBS, by = "row.names")

# Check which degs are significant - 13 DEGS
significant_degs = TPM_sub_cyto_cells_merged %>% filter(padj<=0.05) %>% arrange(desc(log2FoldChange))
significant_degs = column_to_rownames(significant_degs, var = "Row.names")

# Append dataframne to show DEG is significant
TPM_sub_cyto_cells_merged$significance = " "
TPM_sub_cyto_cells_merged = column_to_rownames(TPM_sub_cyto_cells_merged, var = "Row.names")
TPM_sub_cyto_cells_merged$significance <- ifelse(rownames(TPM_sub_cyto_cells_merged) %in% rownames(significant_degs), "*", " ")

# Change rownames to first column
TPM_sub_cyto_cells_merged <- rownames_to_column(TPM_sub_cyto_cells_merged, "genes") # Apply rownames_to_column

# Make a new column combining gene name and significance 
TPM_sub_cyto_cells_merged$row_names_combined <- paste0(TPM_sub_cyto_cells_merged$genes, " ", TPM_sub_cyto_cells_merged$significance)

# Set the new column as row names and remove unnecessary columns
rownames(TPM_sub_cyto_cells_merged) <- TPM_sub_cyto_cells_merged$row_names_combined
TPM_sub_cyto_cells_merged <- TPM_sub_cyto_cells_merged[, -c(1, ncol(TPM_sub_cyto_cells_merged))]
TPM_sub_cyto_cells_merged_heatmap <- TPM_sub_cyto_cells_merged[, 1:13]

# Manually Scale the Data 
TPM_sub_cyto_cells_merged_heatmap.scaled <- pheatmap.scale(TPM_sub_cyto_cells_merged_heatmap)
                                  
# Export the data 
write.csv(TPM_sub_cyto_cells_merged_heatmap.scaled, "/path_to/TPM_sub_cytotoxic_T_cell_merged_heatmap_scaled.csv")

#  Gamma_delta Cell
gamma_delta_cells_cells_cell_expression = rownames(DEGs_AbC_vs_PBS)[rownames(DEGs_AbC_vs_PBS) %in% gamma_delta_cells]

# Subset the TPM data to include only GOI 
TPM_sub_gd_cells = TPM_sub[gamma_delta_cells_cells_cell_expression, ]
TPM_sub_gd_cells = TPM_sub_gd_cells %>% filter(rowSums(across(where(is.numeric)))!=0)

# Determine heatmap annotation
rat_info_sub %>% group_by(Treatment_short)
sample_annotation = rat_info_sub[("Treatment_short")]
rownames(sample_annotation) <- colnames(TPM_sub_gd_cells)
TPM_sub_gd_cells = TPM_sub_gd_cells[do.call(order, c(TPM_sub_gd_cells, list(decreasing=TRUE))),]

# Order the TPM/DEG selected data from high to low
TPM_sub_gd_cells = TPM_sub_gd_cells[do.call(order, c(TPM_sub_gd_cells, list(decreasing=TRUE))),]

# Merge the DEGS/TPM data with DEGs metadata
TPM_sub_gd_cells_merged <- merge(TPM_sub_gd_cells, DEGs_AbC_vs_PBS, by = "row.names")

# Check which degs are significant
significant_degs = TPM_sub_gd_cells_merged %>% filter(padj<=0.05) %>% arrange(desc(log2FoldChange))
significant_degs = column_to_rownames(significant_degs, var = "Row.names")

# Append dataframne to show DEG is significant
TPM_sub_gd_cells_merged$significance = " "
TPM_sub_gd_cells_merged = column_to_rownames(TPM_sub_gd_cells_merged, var = "Row.names")
TPM_sub_gd_cells_merged$significance <- ifelse(rownames(TPM_sub_gd_cells_merged) %in% rownames(significant_degs), "*", " ")

# Change rownames to first column
TPM_sub_gd_cells_merged <- rownames_to_column(TPM_sub_gd_cells_merged, "genes") # Apply rownames_to_column

# Make a new column combining gene name and significance 
TPM_sub_gd_cells_merged$row_names_combined <- paste0(TPM_sub_gd_cells_merged$genes, " ", TPM_sub_gd_cells_merged$significance)

# Set the new column as row names and remove unnecessary columns
rownames(TPM_sub_gd_cells_merged) <- TPM_sub_gd_cells_merged$row_names_combined
TPM_sub_gd_cells_merged <- TPM_sub_gd_cells_merged[, -c(1, ncol(TPM_sub_gd_cells_merged))]
TPM_sub_gd_cells_merged_heatmap <- TPM_sub_gd_cells_merged[, 1:13]

# Manually Scale the Data 
TPM_sub_gd_cells_merged_heatmap.scaled <- pheatmap.scale(TPM_sub_gd_cells_merged_heatmap)

# Export the data 
write.csv(TPM_sub_gd_cells_merged_heatmap.scaled, "/path_to/TPM_sub_gamma_delta_T_cell_merged_heatmap_scaled.csv")


