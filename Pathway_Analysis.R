## Author: Katie Stoker 

# A script to create heatmaps to compare the pathways of interest in PBS vs treated rats. 
# Substitute each pathway of interest into this code.  

# Set up environment
library("dplyr")
library("pheatmap")
library("tidyr")
library("stringr")
library("tidyverse")

# Define the pheatmap.scale function

pheatmap.scale <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  return((x - m) / s)
}

# Import pathway data 

# Data is filtered to include results only with padj <= 0.05 already
pbs_abc_pathways = read.csv("/path_to/Ab_C_vs_PBS_fgsea_results_Reactome_padj_filtered.csv")

# Take a look at the data (8th column is the leading edge)
pbs_abc_pathways[1:5, 1:7]

# Import DEG information
DEGs_AbC_vs_PBS = read.csv("/path_to/DEGs_AbC_vs_PBS.csv", row.names = 1)

# Import the transcript per million data 
TPM = readRDS("/path_to/TPM_counts_collapsed_multiplicates_human_genes.RDS")
rownames(TPM) = TPM$human_symbol
TPM$human_symbol = NULL
TPM$symbol = NULL
colnames(TPM) = gsub("\\.", "-", colnames(TPM))

# Check the pathway  of interest is in imported data
### ADD YOUR PATHWAY OF INTEREST HERE ###
pbs_abc_pathways$pathway[grep("REACTOME_SIGNALING_BY_PDGF", pbs_abc_pathways$pathway)]

# Define pathway of interest
    ### ADD YOUR PATHWAY OF INTEREST HERE ###
POI = "REACTOME_SIGNALING_BY_PDGF"

# Define genes of interest (those in the leading edge)
GOI = pbs_abc_pathways %>% filter(pathway == POI) %>%
  select(leadingEdge)

# Function to split the leadingEdge column and unnest the result 
split_and_unnest <- function(GOI, leadingEdge) {
  GOI %>%
    mutate(gene = strsplit(leadingEdge, " ")) %>%
    unnest(gene)
}

# Apply the function to the sample data
result <- split_and_unnest(GOI, "leadingEdge")
# Remove the comma from the end of each gene 
result$gene <- gsub(",", "", result$gene)
# Display the result
print(result)
GOI = result %>% select(gene)


# Import the rat information (animal ID, treatment etc.)
rat_info = read.csv("/path_to/rat_information.csv")

# Filter for only the PBS and AbC (26) treated rats 
rat_info_sub = rat_info %>% filter(Treatment_short %in% c("PBS", "Ab_C"))
rownames(rat_info_sub) = rat_info_sub$Animal_ID

# State the annotations you want to use in the heatmap
sample_annotation = rat_info_sub[,c("Animal_ID", "Treatment_short")]
sample_annotation$Animal_ID = NULL

# Merge the gene expression data with the list of GOI from the POI
ge = TPM[,]

# Sort the genes from the POI to match those in the TPM data 
sorted_genes = rownames(ge)[rownames(ge) %in% GOI$gene]

# Match animal ID between the datasets (so the rows of sample annotation with the columns of the TPM data)
sorted_animals = rownames(sample_annotation)[rownames(sample_annotation) %in% colnames(ge)]

# Create heatmap of all genes in pathway based on TPM data 
heatmap_table = ge[sorted_genes, sorted_animals,]

# Order in descending order
heatmap_table_sorted = heatmap_table[do.call(order, c(heatmap_table, list(decreasing=TRUE))),]


# Generate heatmaps based on DEGs for each pathway 

# Sort by fold change (descending)
DEGs_AbC_vs_PBS_ordered = DEGs_AbC_vs_PBS %>% arrange(desc(log2FoldChange))

# Find the DEGS which are in the pathway list of genes 
DEGS_in_GOI = rownames(DEGs_AbC_vs_PBS_ordered)[rownames(DEGs_AbC_vs_PBS_ordered) %in% GOI$gene]
DEGS_in_GOI = as.matrix(DEGS_in_GOI)

# FIlter the list of DEGs to include only those from pathway 
DEGs_AbC_vs_PBS_ordered_new = DEGs_AbC_vs_PBS_ordered %>% filter(row.names(DEGs_AbC_vs_PBS_ordered) %in% c(DEGS_in_GOI))

# Select top 50 DEGs for heatmap
DEGs_AbC_vs_PBS_top_50 =  top_n(DEGs_AbC_vs_PBS_ordered_new, n = 50, wt = log2FoldChange)

# Change the rownames (genes) to the first column so datasets can be matched
DEGs_AbC_vs_PBS_top_50 <- rownames_to_column(DEGs_AbC_vs_PBS_top_50, "genes")

# Select TPM data with only the rats of interest
TPM_subset = TPM[, rat_info_sub$Animal_ID]

# Change the rownames (genes) to the first column so datasets can match
TPM_subset <- rownames_to_column(TPM_subset, "genes") # Apply rownames_to_column

# Select the TPM data based on only the top 50 DEGs
TPM_deg_df = left_join(TPM_subset, DEGs_AbC_vs_PBS_top_50, by = "genes")

# Select the top 50 rows based on FC
TPM_sub_50 =  top_n(TPM_deg_df, n = 50, wt = log2FoldChange)
TPM_sub_50 = TPM_sub_50 %>% arrange(desc(log2FoldChange))

# Change the first column back to rownames 
TPM_sub_50 = column_to_rownames(TPM_sub_50, "genes") # Apply rownames_to_column

sorted_animals_degs = rownames(sample_annotation)[rownames(sample_annotation) %in% colnames(TPM_sub_50)]

# Check which degs are significant 
TPM_sub_50_sig = TPM_sub_50 %>% filter(padj<=0.05) %>% arrange(desc(log2FoldChange))

# Append dataframne to show DEG is significant
TPM_sub_50$significance = " "
TPM_sub_50$significance <- ifelse(rownames(TPM_sub_50) %in% rownames(TPM_sub_50_sig), "*", " ")

# Change rownames to first column
TPM_sub_50 <- rownames_to_column(TPM_sub_50, "genes") # Apply rownames_to_column

# Make a new column combining gene name and significance 
TPM_sub_50$row_names_combined <- paste0(TPM_sub_50$genes, " ", TPM_sub_50$significance)

# Set the new column as row names and remove unnecessary columns
rownames(TPM_sub_50) <- TPM_sub_50$row_names_combined
TPM_sub_50 <- TPM_sub_50[, -c(1, ncol(TPM_sub_50))]

# Create heatmap table 
heatmap_table_degs = TPM_sub_50[sorted_animals_degs]

# Manually Scale the Data 
heatmap_table_degs.scaled <- pheatmap.scale(heatmap_table_degs)

# Export the data 
write.csv(heatmap_table_degs.scaled, paste("/path_to/pathway_heatmaps_new/", gsub(" ", "_", POI), "_scaled.csv", sep=""))
