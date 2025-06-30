
################################################################# Packages ############################################
#######################################################################################################################
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(knitr)
library(kableExtra)
library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(SeuratWrappers)
library(DoubletFinder)
library(SoupX)


###################################################### Reading Seurat Clustering object#############################################
####################################################################################################################################
rm(list=ls())
setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/2_Seurat/Integrate/Combined/Separate/comparison-c9-Vs-all-find-markers/integrated_second-clustering_excluding-c9/")
sample<- readRDS("clustering-r5-pc21_SCT.rds")
dim(sample)
DefaultAssay(sample)<-"RNA"
#23756 25765
DimPlot(sample, reduction = "umap", label = TRUE, label.size = 8, pt.size = 1) 
DimPlot(sample, reduction = "umap", label = TRUE, label.size = 8, pt.size = 1, split.by = "orig.ident") 



####################################### UMI/Gene distribution ####################################################
##################################################################################################################

# Generate a violin plot
vln_plot <- VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Save the plot as a PDF
pdf("/depot/tlratlif/data/LSC_scRNAseq_2023/Contaminant-suspension/violin_plot.pdf", width = 8, height = 6)  # Specify width and height
print(vln_plot)
dev.off()  # Close the PDF device

# Generate a scatter plot
scatter_plot <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Save the plot as a PDF
pdf("/depot/tlratlif/data/LSC_scRNAseq_2023/Contaminant-suspension/scatter_plot.pdf", width = 8, height = 6)
print(scatter_plot)
dev.off()

# all together 

v2<-VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol=4)

# Save the plot as a PDF
pdf("/depot/tlratlif/data/LSC_scRNAseq_2023/Contaminant-suspension/violin_plot_all_4.pdf", width = 8, height = 6)  # Specify width and height
print(v2)
dev.off()  # Close the PDF device


##################################### Mitochondrial content ######################################################
##################################################################################################################


# Generate vlnplot showing mitco content
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
mit<-VlnPlot(sample, features = "percent.mt")

pdf("/depot/tlratlif/data/LSC_scRNAseq_2023/Contaminant-suspension/mitochondrial_plot.pdf", width = 8, height = 6)
print(mit)
dev.off()


##################################### Doublet identification ######################################################
##################################################################################################################


library(Seurat)
library(DoubletFinder)

sample <- NormalizeData(sample)
sample <- FindVariableFeatures(sample)
sample <- ScaleData(sample)
sample <- RunPCA(sample)

# Estimate expected number of doublets
sweep.res.list <- paramSweep_v3(sample, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Choose the best pK value (highest bcmvn peak)
best.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Run DoubletFinder with expected doublet rate (adjust rate based on your dataset)
nExp <- round(0.075 * nrow(sample))  # Example: 7.5% doublet rate
seurat_object <- doubletFinder_v3(sample, PCs = 1:10, pN = 0.25, pK = best.pK, nExp = nExp, reuse.pANN = FALSE, sct = FALSE)

##  Stacked bar plot



# Make sure the identity is set to clusters
Idents(seurat_object) <- "seurat_clusters"

# Find the name of the classification column added by DoubletFinder
df_col <- grep("DF.classifications", colnames(seurat_object@meta.data), value = TRUE)

# Create a summary table: counts of Singlet/Doublet per cluster
df_summary <- seurat_object@meta.data %>%
  dplyr::select(seurat_clusters, !!df_col) %>%
  group_by(seurat_clusters, !!sym(df_col)) %>%
  summarise(count = n(), .groups = "drop")

# Plot the bar graph
ggplot(df_summary, aes(x = seurat_clusters, y = count, fill = !!sym(df_col))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Doublets and Singlets Across Clusters",
       x = "Cluster",
       y = "Cell Count",
       fill = "Classification") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))


################################# Cell cycle score ###############################################################
##################################################################################################################
# Load Seurat's built-in cell cycle gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


s.genes <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG")
g2m.genes <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TOP2A", "TPX2", "MKI67")


# Perform cell cycle scoring
seurat_object <- CellCycleScoring(
    object = sample,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE # Optionally set cell cycle phase as identity class
)



# visualize cell cycle scores

cell.cyle<-VlnPlot(seurat_object, features = c("S.Score", "G2M.Score"), group.by = "seurat_clusters")
pdf("/depot/tlratlif/data/LSC_scRNAseq_2023/Contaminant-suspension/cell_cycle.pdf", width = 8, height = 6)
print(cell.cyle)
dev.off()


# Calculate the proportion of cells in each cell cycle phase across clusters:
table1<-table(seurat_object$Phase, seurat_object$seurat_clusters) %>%
  prop.table(margin = 2) %>%
  as.data.frame()
  
setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/Contaminant-suspension/")
write.table(table1,file = "cell.scores.txt", sep = "\t", quote = F, row.names = T)


