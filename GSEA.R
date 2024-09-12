## Author: Roman Laddach
# A script to perform functional analysis on the list of DEGs. Genes are for human after conversion from rat.

# Set up environment 
library("dplyr")
library("ggplot2")
library("fgsea") 
library("gprofiler2") 

# Load the DEGs
DEGs = read.csv("/path_to/DEGs_AbC_vs_PBS.csv", row.names = 1)

# Preprocessing  

# Remove rows with missing log2FoldChange  
sum(is.na(DEGs$log2FoldChange))
DEGs = DEGs[!is.na(DEGs$log2FoldChange),]
sum(is.na(DEGs$log2FoldChange)) # Sanity check

# GSEA  

# Order the genes 
de = data.frame(DEGs)
de = de[order(de$log2FoldChange, decreasing = T),]

# Extract the list of gene names with their respective foldchange
ranks = setNames(de$log2FoldChange,rownames(de))



## Hallmark ##
hallmark = gmtPathways("/path_to/h.all.v2023.1.Hs.symbols.gmt") # hallmark
set.seed(13)
hallmark_fgseaRes = fgseaMultilevel(hallmark, ranks, nPermSimple = 10000, eps = 0)

# Plot the results 
ggplot(hallmark_fgseaRes, aes(y=reorder(hallmark, NES), x=NES, fill=padj)) +geom_col() 

# Filter for signifance and plot
hallmark_fgseaRes %>% filter(padj <=0.05) %>% ggplot(aes(y=reorder(hallmark, NES), x=NES, fill=padj)) +geom_col() 



## Reactome ##
reactome = gmtPathways("/path_to/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
set.seed(13)
reactome_fgseaRes = fgseaMultilevel(reactome, ranks, nPermSimple = 10000, eps = 0)

# Filter for signifance and plot
reactome_fgseaRes %>% filter(padj <=0.05) %>% ggplot(aes(y=reorder(pathway, NES), x=NES, fill=padj)) +geom_col() 
