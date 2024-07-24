# GENE SET ENRICHMENT ANALYSIS -------------------------------------------------
# NAME: G. Brett Moreau
# DATE: July 9, 2024

# PACKAGES ---------------------------------------------------------------------
#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # Using version 1.3.2

#install.packages("msigdbr")
library(msigdbr)
packageVersion("msigdbr") # Using version 7.5.1

#BiocManager::install("fgsea")
library(fgsea)
packageVersion("fgsea") # Using version 1.20.0

#BiocManager::install("DESeq2")
library(DESeq2)
packageVersion("DESeq2") # Using version 1.34.0

#BiocManager::install("plyranges")
library(plyranges)
packageVersion("plyranges") # Using version 1.14.0

#install.packages("pheatmap")
library(pheatmap)
packageVersion("pheatmap") # Using version 1.0.12

#install.packages("ggrepel")
library(ggrepel)
packageVersion("ggrepel") # Using version 0.9.1

# INTRODUCTION -----------------------------------------------------------------

# Differentially expressed genes were identified between Alum and Liposome 
# adjuvent groups. I will now use the DEG analysis results to perform Gene Set
# Enrichment Analysis (GSEA), identifying biological pathways enriched in either
# group.




# IMPORT AND ORGANIZATION OF RNA SEQUENCING RESULTS ----------------------------

# Import DESeq2 results files
load("../results/DESeq2 results_adjuvent.RData")

# Import gene name information for appending results files.
gene_names <- read_gff("../data/Mus_musculus.GRCm39.109.gtf.gz") 
gene_names <- plyranges::select(gene_names, gene_id, gene_name)
gene_names <- as.data.frame(gene_names)
gene_names <- gene_names %>%
  dplyr::select(gene_name, gene_id) %>%
  dplyr::distinct()


# Add gene names to results files
res_adjuvent_table <- as.data.frame(res_adjuvent)
res_adjuvent_table$gene_id <- rownames(res_adjuvent)
res_adjuvent_table <- res_adjuvent_table %>%
  dplyr::left_join(., gene_names, by = "gene_id") %>%
  dplyr::relocate(gene_id, gene_name, .before = baseMean)




# ORGANIZE RANKED GENE LIST FOR GSEA -------------------------------------------

# Check for missing Wald statistic values
table(is.na(res_adjuvent_table$stat)) # no missing values

# Check for duplicate gene names
table(duplicated(res_adjuvent_table$stat)) # 18 duplicates

# ENSEMBL IDs with identical names will be averaged together.

# Generate ranked gene list.
ranked_adjuvent <- res_adjuvent_table %>%
  dplyr::select(gene_name, stat) %>%
  na.omit() %>%
  group_by(gene_name) %>%
  summarize(stat = mean(stat)) %>%
  deframe()

head(ranked_adjuvent, 20) # Gene list with Wald statistic for each gene.




# ORGANIZATION OF GENE SETS ----------------------------------------------------

# GSEA will be performed with the fgsea package. I will look for enrichment 
# using the Hallmark, GO: Biological Processes, and KEGG gene sets from the 
# MSigDB Collection. 


### HALLMARK ### 

# Load Hallmark gene sets for mouse
h_gene_set <- msigdbr(species = "Mus musculus", category = "H")

# Organize Hallmark pathway set for GSEA
h_gene_set_list <- split(x = h_gene_set$gene_symbol, 
                         f = h_gene_set$gs_name)

View(h_gene_set_list) # List with genes for each gene set name.



### GENE ONTOLOGY: BIOLOGICAL PROCESSES ###

# Load GO:BP gene sets
GO_BP_set <- msigdbr(species = "Mus musculus", 
                     category = "C5", 
                     subcategory = "GO:BP")

GO_BP_set_list <- split(x = GO_BP_set$gene_symbol, 
                        f = GO_BP_set$gs_name)


### KEGG ###

# Load KEGG gene sets
KEGG_set <- msigdbr(species = "Mus musculus", 
                    category = "C2", 
                    subcategory = "KEGG")

KEGG_set_list <- split(x = KEGG_set$gene_symbol, 
                       f = KEGG_set$gs_name)



# GSEA: HALLMARK GENE SET ------------------------------------------------------

# Run GSEA
fgsea_adjuvent_Hallmark <- fgsea(pathways = h_gene_set_list, 
                                 stats = ranked_adjuvent, 
                                 eps = 0)


# Collapse leadingEdge column from list for table printing
fgsea_adjuvent_Hallmark_table <- fgsea_adjuvent_Hallmark %>%  
  rowwise() %>%
  mutate(leadingEdge = paste(leadingEdge, collapse=',')) %>%
  ungroup()

# Count number of commas in leading edge column as proxy for number of genes.
fgsea_adjuvent_Hallmark_table$enriched_genes <- 
  lengths(gregexpr(",", fgsea_adjuvent_Hallmark_table$leadingEdge)) + 1

# Calculate the percent of pathway genes enriched for each pathway
fgsea_adjuvent_Hallmark_table$percent_coverage <- 
  fgsea_adjuvent_Hallmark_table$enriched_genes / 
  fgsea_adjuvent_Hallmark_table$size * 100



### OVERVIEW OF SIGNIFICANTLY ENRICHED PATHWAYS ###

# Count number of DEGs overall
sum(fgsea_adjuvent_Hallmark_table$padj < 0.05, na.rm = TRUE) # 26 gene sets

# Count gene sets per group
fgsea_adjuvent_Hallmark_table %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(NES > 0) %>%
  dplyr::select(pathway) # 16 gene sets enriched in Liposome condition.


fgsea_adjuvent_Hallmark_table %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(NES < 0) %>% 
  dplyr::select(pathway) # 10 gene sets enriched in Alum condition



# Print results files (remove comments to save tables)
#write.csv(fgsea_adjuvent_Hallmark_table, 
#          file = "../results/tables/GSEA adjuvent_Hallmark.csv", 
#          row.names = FALSE)




# GSEA: GO BP GENE SET ------------------------------------------------------------

# Run GSEA
fgsea_adjuvent_GO_BP <- fgsea(pathways = GO_BP_set_list, 
                                 stats = ranked_adjuvent, 
                                 eps = 0)


# Collapse leadingEdge column from list for table printing
fgsea_adjuvent_GO_BP_table <- fgsea_adjuvent_GO_BP %>%  
  rowwise() %>%
  mutate(leadingEdge = paste(leadingEdge, collapse=',')) %>%
  ungroup()

# Count number of commas in leading edge column as proxy for number of genes.
fgsea_adjuvent_GO_BP_table$enriched_genes <- 
  lengths(gregexpr(",", fgsea_adjuvent_GO_BP_table$leadingEdge)) + 1

# Calculate the percent of pathway genes enriched for each pathway
fgsea_adjuvent_GO_BP_table$percent_coverage <- 
  fgsea_adjuvent_GO_BP_table$enriched_genes / 
  fgsea_adjuvent_GO_BP_table$size * 100



### OVERVIEW OF SIGNIFICANTLY ENRICHED PATHWAYS ###

# Count number of DEGs overall
sum(fgsea_adjuvent_GO_BP_table$padj < 0.05, na.rm = TRUE) # 790 gene sets

# Count gene sets per group
fgsea_adjuvent_GO_BP_table %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(NES > 0) %>%
  dplyr::select(pathway) # 664 gene sets enriched in Liposome condition.


fgsea_adjuvent_GO_BP_table %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(NES < 0) %>% 
  dplyr::select(pathway) # 116 gene sets enriched in Alum condition



# Print results files (remove comments to save tables)
#write.csv(fgsea_adjuvent_GO_BP_table, 
#          file = "../results/tables/GSEA adjuvent_GO_BP.csv", 
#          row.names = FALSE)




# GSEA: KEGG GENE SET ----------------------------------------------------------

# Run GSEA
fgsea_adjuvent_KEGG <- fgsea(pathways = KEGG_set_list, 
                              stats = ranked_adjuvent, 
                              eps = 0)


# Collapse leadingEdge column from list for table printing
fgsea_adjuvent_KEGG_table <- fgsea_adjuvent_KEGG %>%  
  rowwise() %>%
  mutate(leadingEdge = paste(leadingEdge, collapse=',')) %>%
  ungroup()

# Count number of commas in leading edge column as proxy for number of genes.
fgsea_adjuvent_KEGG_table$enriched_genes <- 
  lengths(gregexpr(",", fgsea_adjuvent_KEGG_table$leadingEdge)) + 1

# Calculate the percent of pathway genes enriched for each pathway
fgsea_adjuvent_KEGG_table$percent_coverage <- 
  fgsea_adjuvent_KEGG_table$enriched_genes / 
  fgsea_adjuvent_KEGG_table$size * 100



### OVERVIEW OF SIGNIFICANTLY ENRICHED PATHWAYS ###

# Count number of DEGs overall
sum(fgsea_adjuvent_KEGG_table$padj < 0.05, na.rm = TRUE) # 52 gene sets

# Count gene sets per group
fgsea_adjuvent_KEGG_table %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(NES > 0) %>%
  dplyr::select(pathway) # 37 gene sets enriched in Liposome condition.


fgsea_adjuvent_KEGG_table %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(NES < 0) %>% 
  dplyr::select(pathway) # 15 gene sets enriched in Alum condition



# Print results files (remove comments to save tables)
#write.csv(fgsea_adjuvent_KEGG_table, 
#          file = "../results/tables/GSEA adjuvent_KEGG.csv", 
#          row.names = FALSE)

