# IMPORT AND ANNOTATION OF KALLISTO DATA ---------------------------------------
# NAME: G. Brett Moreau
# DATE: July 9, 2024

# PACKAGES ---------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
packageVersion("BiocManager") # I'm using version 1.30.22

#install.packages("readxl")
library(readxl)
packageVersion("readxl") # Using version 1.4.3

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # Using version 2.0.0

#BiocManager::install("rhdf5")
library(rhdf5)
packageVersion("rhdf5") # Using version 2.46.1

#BiocManager::install("tximport")
library(tximport)
packageVersion("tximport") # Using version 1.30.0

#BiocManager::install("ensembldb")
library(ensembldb)
packageVersion("ensembldb") # Using version 2.26.0

#BiocManager::install("DESeq2")
library(DESeq2)
packageVersion("DESeq2") # Using version 1.42.0

#BiocManager::install("plyranges")
library(plyranges)
packageVersion("plyranges") # Using version 1.22.0

#install.packages("RColorBrewer")
library(RColorBrewer)
packageVersion("RColorBrewer") # Using version 1.1.3

#install.packages("pheatmap")
library(pheatmap)
packageVersion("pheatmap") # Using version 1.0.12

# INTRODUCTION -----------------------------------------------------------------

# Bulk RNA sequencing was performed by Novogene on RNA isolated from mice given
# different vaccine regimens, then challenged with C. difficile. Read quality 
# was analyzed using FastQC, Illumina adapters were removed using BBDuk, and 
# reads were pseudoaligned with the murine genome using Kallisto. Mapped reads 
# will now be imported into R and analyzed using DESeq2.




# ORGANIZE METADATA ------------------------------------------------------------

# Import metadata and format for DESeq2
metadata <- read.csv("../data/metadata.csv")

metadata <- tibble::column_to_rownames(metadata, "SampleID")
metadata <- cbind("SampleID" = rownames(metadata), metadata)

metadata <- metadata %>%
  dplyr::mutate(Group = factor(Group, levels = c("Control", "Alum", "Liposome"))) %>%
  dplyr::mutate(Cage_Number = factor(Cage_Number)) %>%
  dplyr::mutate(SampleID = factor(SampleID))




# IMPORT KALLISTO DATA INTO R  -------------------------------------------------

# Generate file path for Kallisto files
path <- file.path("../results/kallisto", 
                  metadata$SampleID, 
                  "abundance.tsv")

file.exists(path) # All files present
names(path) <- metadata$SampleID


# Now I need to make a tx2gene file, which will allow us to match transcript and 
# gene identifiers during import of Kallisto abundance information.Because I 
# made my reference index from the Mus musculis cDNA FASTA on the Ensembl 
# website, I need to use a matching transcript database for matching gene 
# identifiers during import. I will therefore use the GTF from the Ensemble FTP 
# site (VERSION 109, MATCHING THE INDEX VERSION) for generating the transcript 
# database.

# Make Tx2gene file for Tximport
TxDB<- GenomicFeatures::makeTxDbFromGFF("../data/Mus_musculus.GRCm39.109.gtf.gz",
                                        format = "gtf")
k <- keys(TxDB, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(TxDB, k, "GENEID", "TXNAME")


# Import Kallisto data. Will also estimate counts using length-adjusted 
# abundance estimates.
txi_kallisto <- tximport(path, 
                         type = "kallisto", 
                         tx2gene = tx2gene, 
                         countsFromAbundance = "lengthScaledTPM",
                         ignoreTxVersion = TRUE)

View(txi_kallisto$counts) # Scaled counts with gene identifiers as row names.

# There are two methods of data normalization that are recommended by the makers
# of tximport before handing off to differential gene expression pipelines. The
# method used here scales counts using the average transcript length over 
# samples and the library size, then hands these adjusted counts (in the 
# 'counts' table) directly to the DGE pipeline of choice. I'll use DESeq2 for
# differential gene expression analysis.




# HAND OFF COUNTS TO DESEQ2 ----------------------------------------------------

# Confirm identical names for metadata and counts.
all(rownames(metadata) == colnames(txi_kallisto$counts)) # Identical names

# Make DESeq2 data set
dds <- DESeqDataSetFromTximport(txi = txi_kallisto, 
                                colData = metadata, 
                                design = ~ Group)

# Filter low abundance genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Generate DESeq object for differential gene expression
dds <- DESeq(dds)




# EVALUATE GROUP DIFFERENCES BY PRINCIPAL COMPONENT ANALYSIS -------------------

# Log transform and organize count data for PCA
rld <- rlog(dds, blind = TRUE)
plotpca <- plotPCA(rld, intgroup = c("Group"), returnData = TRUE)
percentVar <- round(100 * attr(plotpca, "percentVar"))

# Plot PCA
ggplot(data = plotpca, aes(x = PC1, y = PC2,
                           #color = metadata$Group,
                           label = metadata$SampleID
)) + 
  geom_point(aes(color = Group), size = 6) +
  geom_text(hjust=0, vjust=0) +
  labs(color = "Group") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        aspect.ratio = 1)

#ggsave("../results/figures/PCA plot_labeled.png", 
#      width = 7, dpi = 300)

# Control (unvaccinated) samples cluster distinctly from vaccinated samples, as
# expected.


### HEATMAP OF WITHIN-SAMPLE VARIABILITY ###

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$SampleID, " (", rld$Group, ")")
colnames(sampleDistMatrix) <- paste(rld$SampleID, " (", rld$Group, ")")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#png("../results/figures/sample distances.png", 
#    width = 750, height = 750)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
#dev.off()




# REMOVAL OF UNVACCINATED SAMPLES AND REMAKE OF DESEQ2 OBJECT ------------------

# I will remove the unvaccinated control samples and look at separation between
# adjuvent groups.

# Remove sample from metadata
metadata_no_controls <- metadata %>%
  dplyr::filter(Group != "Control") %>%
  dplyr::mutate(Group = recode(Group, "Alum" = "TcdB-Alum", "Liposome" = "TcdB-LS")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("TcdB-Alum", "TcdB-LS"))) %>%
  dplyr::mutate(SampleID = factor(SampleID))


# Remove control samples from path list.
control_paths <- c("../results/kallisto/M1/abundance.tsv", 
                   "../results/kallisto/M2/abundance.tsv",
                   "../results/kallisto/M3/abundance.tsv",
                   "../results/kallisto/M4/abundance.tsv")
?case_match()

# Import non-control samples
path_no_controls <- path[!path %in% control_paths]

txi_kallisto_no_controls <- tximport(path_no_controls, 
                                     type = "kallisto", 
                                     tx2gene = tx2gene, 
                                     countsFromAbundance = "lengthScaledTPM", 
                                     ignoreTxVersion = TRUE)


# Confirm identical names for metadata and counts.
all(rownames(metadata_no_controls) == colnames(txi_kallisto_no_controls$counts)) 


# Make DESeq2 data set
dds_no_controls <- DESeqDataSetFromTximport(txi = txi_kallisto_no_controls, 
                                            colData = metadata_no_controls, 
                                            design = ~ Group)


# Filter low abundance genes
keep_no_controls <- rowSums(counts(dds_no_controls)) >= 10
dds_no_controls <- dds_no_controls[keep_no_controls,]

# Generate DESeq object for differential gene expression
dds_no_controls <- DESeq(dds_no_controls)




# VISUALIZATION OF SAMPLE CLUSTERING: VACCINE SAMPLES ONLY ---------------------


# Log transform and organize count data for PCA
rld_no_controls <- rlog(dds_no_controls, blind = TRUE)
plotpca_no_controls <- plotPCA(rld_no_controls, intgroup = c("Group"), returnData = TRUE)
percentVar_no_controls <- round(100 * attr(plotpca_no_controls, "percentVar"))

# Plot PCA
ggplot(data = plotpca_no_controls, aes(x = PC1, y = PC2,
                           label = metadata_no_controls$SampleID
)) + 
  geom_point(aes(color = Group), size = 6) +
  #geom_text(hjust=0, vjust=0) +
  labs(color = "Group") +
  xlab(paste0("PC1: ",percentVar_no_controls[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_no_controls[2],"% variance")) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        aspect.ratio = 1)

#ggsave("../results/figures/PCA plot_noControls.png", 
#       width = 7, dpi = 300)


### HEATMAP OF WITHIN-SAMPLE VARIABILITY ###

sampleDists_no_controls <- dist(t(assay(rld_no_controls)))

sampleDistMatrix_no_controls <- as.matrix(sampleDists_no_controls)
rownames(sampleDistMatrix_no_controls) <- paste(rld_no_controls$SampleID, " (", rld_no_controls$Group, ")")
colnames(sampleDistMatrix_no_controls) <- paste(rld_no_controls$SampleID, " (", rld_no_controls$Group, ")")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#png("../results/figures/sample distances_noControls.png", 
#    width = 750, height = 750)
pheatmap(sampleDistMatrix_no_controls,
         clustering_distance_rows = sampleDists_no_controls,
         clustering_distance_cols = sampleDists_no_controls,
         col = colors)
#dev.off()




# GENERATE RESULTS FILES FOR COMPARISONS OF INTEREST ---------------------------

# I will focus on the Adjuvent comparison (Alum vs Liposome) for now.

# Generate results files for each comparison.
res_adjuvent <- results(dds_no_controls, contrast = c("Group", 
                                            "TcdB-LS", 
                                            "TcdB-Alum"))




# QUALITY CONTROL FOR RNA SEQUENCING DATA --------------------------------------

### DISPERSION ANALYSIS ###

# Plot Dispersion Estimates (remove comments to save plots)
#png("../results/figures/dispersion estimates_noControls.png", 
#    width = 750, height = 750)
plotDispEsts(dds_no_controls)
#dev.off()

# The final dispersion estimates are shrunk compared to the gene-wise estimates.
# This is the expected result from this analysis.


### MA-PLOT ###
# Generate and save MA-plots (remove comments to save plots)

#png("../results/figures/MA_adjuvent.png", 
#    width = 750, height = 750)
DESeq2::plotMA(res_adjuvent)
#dev.off()


# Save only results files for downstream analysis
rm(list = setdiff(ls(), c("res_adjuvent", 
                          "rld_no_controls", 
                          "dds_no_controls",
                          "metadata_no_controls")))

#save.image(file = "../results/DESeq2 results_adjuvent.RData")



