---
title: "Post-hoc"
author: "Nikhila"
date: "`r Sys.Date()`"
output: html_document
---


```{r library, warning=FALSE, message=FALSE, echo=FALSE, tidy=TRUE, eval=TRUE}

library(Seurat)
library(Seurat)
library(dplyr)
library(tidyverse)
library(knitr)
library(kableExtra)
library(ggthemes)
library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(glmGamPoi)
library(clustree)
library(ggraph)
library(SingleR)
library(xlsx)
library(pheatmap)
library(dplyr)
library(AUCell)
library(GSEABase)

```



```{r load_data, warning=FALSE, message=FALSE, echo=FALSE, tidy=TRUE, eval = T}

setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/2_Seurat/Integrate/Combined/Separate")


name = "clustering-r6-naiv-infl"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )
DefaultAssay(sample) <- "RNA"

title = paste0("Input sample having 14 clusters: nPCs = ", 30,"resolution = ", 0.6) 


pdf("initial_clustering_naiv_infl_umap.pdf", width = 10, height = 10)
DimPlot(sample , reduction = "umap", label = TRUE, label.size = 6, split.by = "expt.type")
dev.off()
sample$expt.type <- factor(sample$expt.type, levels = c("Naive", "Inflamed"))

```



```{r aucell, warning=FALSE, message=FALSE, echo=FALSE, tidy=TRUE, eval = T}


setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/2_Seurat/Integrate/Combined/Separate/aucell")



# 1. Extract expression matrix from Seurat object
exprMatrix <- as.matrix(GetAssayData(sample, slot = "data"))  # log-normalized data

# 2. Build gene rankings per cell
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = TRUE, nCores = 1)

# 3. Define your gene sets (use your own marker sets)
geneSets <- list(
  Basal_response = c("Krt5","Krt14", "Trp63"),         # example basal markers
  Luminal_Signature = c("Krt8", "Krt18", "Cd24a")       # example luminal markers
)

# 4. Calculate AUC for each cell and gene set
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# 5. Inspect available AUC data
rownames(getAUC(cells_AUC))  # should show: "Basal_response", "Luminal_Signature"

# 6. Transfer AUC scores to Seurat meta.data
auc_matrix <- as.data.frame(t(getAUC(cells_AUC)))
common_cells <- intersect(rownames(sample@meta.data), rownames(auc_matrix))
sample@meta.data[common_cells, "Basal_Response"] <- auc_matrix[common_cells, "Basal_response"]
sample@meta.data[common_cells, "Luminal_Signature"] <- auc_matrix[common_cells, "Luminal_Signature"]



# Reverse the factor levels to swap default colors
sample$expt.type <- factor(sample$expt.type, levels = c("Inflamed", "Naive"))

# Plot with default colors swapped
pdf("basal_au-cell.pdf", width = 10, height = 10)
VlnPlot(
  sample,
  features = "Basal_Response",
  group.by = "seurat_clusters",
  split.by = "expt.type",
  pt.size = 0.1
) + scale_fill_manual(values = c(Naive = "#00BFC4", Inflamed = "#F8766D"))


dev.off()

pdf("luminal-au-cell.pdf", width = 10, height = 10)
VlnPlot(
  sample,
  features = "Luminal_Signature",
  group.by = "seurat_clusters",
  split.by = "expt.type",
  pt.size = 0.1
) + scale_fill_manual(values = c(Naive = "#00BFC4", Inflamed = "#F8766D"))
dev.off()


```




```{r qc-stack, warning=FALSE, message=FALSE, echo=FALSE, tidy=TRUE, eval = T}

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# Extract relevant metadata
qc_data <- sample@meta.data[, c("expt.type", "nFeature_RNA", "nCount_RNA", "percent.mt")]
qc_data$CellID <- rownames(qc_data)

# Bin percent.mt
qc_data$mt_bin <- cut(
  qc_data$percent.mt,
  breaks = c(0, 5, 10, 20, 100),
  labels = c("0-5%", "5-10%", "10-20%", ">20%"),
  include.lowest = TRUE
)

# Bin nFeature_RNA
qc_data$feature_bin <- cut(
  qc_data$nFeature_RNA,
  breaks = c(0, 500, 1000, 2000, 5000, 10000),
  labels = c("<500", "500-1000", "1000-2000", "2000-5000", ">5000"),
  include.lowest = TRUE
)

# Bin nCount_RNA
qc_data$count_bin <- cut(
  qc_data$nCount_RNA,
  breaks = c(0, 1000, 5000, 10000, 25000, 50000, Inf),
  labels = c("<1K", "1K-5K", "5K-10K", "10K-25K", "25K-50K", ">50K"),
  include.lowest = TRUE
)

# Summarize cell counts per bin per condition
mt_summary <- qc_data %>%
  group_by(expt.type, mt_bin) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(metric = "percent.mt", bin = mt_bin)

feature_summary <- qc_data %>%
  group_by(expt.type, feature_bin) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(metric = "nFeature_RNA", bin = feature_bin)

count_summary <- qc_data %>%
  group_by(expt.type, count_bin) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(metric = "nCount_RNA", bin = count_bin)

# Combine all summaries
qc_summary <- bind_rows(mt_summary, feature_summary, count_summary)

# Plot: Faceted stacked bar plots by metric

setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/2_Seurat/Integrate/Combined/Separate/aucell")
pdf("qc-stacked-bar-plots.pdf", width = 10, height = 10)

ggplot(qc_summary, aes(x = expt.type, y = cell_count, fill = bin)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "QC Metrics Distribution by Condition",
    x = "Condition (expt.type)",
    y = "Number of Cells",
    fill = "Value Bin"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()

```






# Load required libraries
library(Seurat)
library(UCell)

# Define your gene signatures
gene_signatures <- list(
  Basal = c("Krt5", "Krt14", "Trp63"),        # Replace with real basal markers
  Luminal = c("Krt8", "Krt18", "Cd24a")     # Replace with real luminal markers
)

# Run UCell scoring and add to metadata
# Correct function call without named 'object' argument
sample <- AddModuleScore_UCell(
  sample,                      # Seurat object passed as first unnamed argument
  features = gene_signatures,
  assay = "RNA"                # Replace with correct assay if needed
)


# Check added columns
head(sample@meta.data[, c("Basal", "Luminal")])

setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/2_Seurat/Integrate/Combined/Separate/aucell")
pdf("luminal-module-scores.pdf", width = 10, height = 10)
VlnPlot(
  sample,
  features = "Luminal_Signature",
  group.by = "seurat_clusters",
  split.by = "expt.type",
  pt.size = 0.1
) + scale_fill_manual(values = c(Naive = "#00BFC4", Inflamed = "#F8766D"))


dev.off()

setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/2_Seurat/Integrate/Combined/Separate/aucell")
pdf("basal-module-scores.pdf", width = 10, height = 10)
VlnPlot(
    sample,
    features = "Basal_Response",
    group.by = "seurat_clusters",
    split.by = "expt.type",
    pt.size = 0.1
) + scale_fill_manual(values = c(Naive = "#00BFC4", Inflamed = "#F8766D"))


dev.off()





