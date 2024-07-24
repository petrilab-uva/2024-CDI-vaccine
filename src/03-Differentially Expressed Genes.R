# DIFFERENTIALLY EXPRESSED GENE ANALYSIS ---------------------------------------
# NAME: G. Brett Moreau
# DATE: July 9, 2024

# PACKAGES ---------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
packageVersion("BiocManager") # I'm using version 1.30.22

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # Using version 2.0.0

#BiocManager::install("DESeq2")
library(DESeq2)
packageVersion("DESeq2") # Using version 1.42.0

#BiocManager::install("plyranges")
library(plyranges)
packageVersion("plyranges") # Using version 1.22.0

#install.packages("pheatmap")
library(pheatmap)
packageVersion("pheatmap") # Using version 1.0.12


# INTRODUCTION -----------------------------------------------------------------

# Mapped reads have been imported into R and counts normalized using DESeq2. I 
# will now investigate differentially expressed genes between different adjuvent
# treatment groups.




# GENERATE DIFFERENTIALLY EXPRESSED GENES TABLE --------------------------------

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


### OVERVIEW OF DIFFERENTIALLY EXPRESSED GENES ###

# Count number of DEGs overall
sum(res_adjuvent_table$padj < 0.05, na.rm = TRUE) # 764 DEGs

# Count Number of DEGs per group
res_adjuvent_table %>%
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>%
  dplyr::select(gene_name) # 435 DEGs upregulated in Liposome condition

res_adjuvent_table %>%
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange < 0) %>%
  dplyr::select(gene_name) # 329 DEGs upregulated in Alum condition


# Print results files (remove comments to save tables)
#write.csv(res_adjuvent_table, 
#          file = "../results/tables/LIPOSOMEvsALUM_DEG_table.csv", 
#          row.names = FALSE)




# VOLCANO PLOT OF DIFFERENTIALLY EXPRESSED GENES -------------------------------

# Generate column for colors
res_adjuvent_volcano <- res_adjuvent_table %>%
  mutate(Significance = ifelse(padj > 0.05, "Not Significant", 
                               ifelse(padj < 0.05 & log2FoldChange >= 0, 
                                      "Upregulated", 
                                      ifelse(padj < 0.05 & log2FoldChange <= 0,
                                             "Downregulated", "No FC")))) %>%
  mutate(Significance = factor(Significance, levels = c("Upregulated",
                                                        "Downregulated", 
                                                        "Not Significant")))


# Plot volcano plot
ggplot(data = res_adjuvent_volcano, aes(x = log2FoldChange,
                                            y = -log10(padj), 
                                            color = Significance)) +
  geom_point(size = 5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  labs(x = "Log2(Fold Change)", y = "Log10(p-value)") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "gray")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title.x = element_text(size = 40), 
        axis.title.y = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))

#ggsave("../results/figures/volcano_adjuvent.png", width = 10, dpi = 300)




# HEATMAPS OF TOP DIFFERENTIALLY EXPRESSED GENES -------------------------------

# Select 50 most differentially expressed genes.
heatmap_adjuvent <- 
  assay(rld_no_controls)[head(order(res_adjuvent$padj), 50), ]


# Organize data set to remove unnecessary samples and add gene names.
heatmap_adjuvent <- as.data.frame(heatmap_adjuvent)
heatmap_adjuvent$gene_id <- rownames(heatmap_adjuvent)
heatmap_adjuvent <- heatmap_adjuvent %>%
  dplyr::left_join(., gene_names, by = "gene_id")


# Convert gene names to row names
rownames(heatmap_adjuvent) <- heatmap_adjuvent$gene_name
heatmap_adjuvent <- dplyr::select(heatmap_adjuvent, -gene_id, -gene_name)


# Convert count data to change in expression across the average of all samples.
heatmap_adjuvent <- heatmap_adjuvent - rowMeans(heatmap_adjuvent)
heatmap_adjuvent <- as.matrix(heatmap_adjuvent)


# Format metadata to add to plot
key <- metadata_no_controls %>% 
  dplyr::select(Group) %>%
  dplyr::mutate(Group = factor(Group, levels = c("TcdB-Alum", "TcdB-LS")))


# Plot heatmap
#png("../results/figures/heatmap_adjuvent_top50.png", 
#    width = 7, height = 10, units = "in", res = 300)
pheatmap(heatmap_adjuvent, 
         annotation_col = key)
#dev.off()

