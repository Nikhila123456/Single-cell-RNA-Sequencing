
#Comparison analysis of multiple datasets with different cell type compositions

library(CellChat)
library(patchwork)

library(Seurat)
library(future)
library(future.apply)

library(ggplot2)
library(tidyr)
library(dplyr)
library(limma)

library(scales)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(grid)

######-------Cellchat object for Stem 


#Reading RDS object

setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/2_Seurat/Integrate/Combined/Separate/aire-comparison/collagen-based-filtering/")
stem = readRDS("Aire-1-final-r1-15-pca_SCT.rds")
DefaultAssay(stem)<-"RNA"
table(stem$seurat_clusters)
#0  1  2 
#66 44 34 


#Adding meta column:cc_labels


cell.type="Stem"
stem <- AddMetaData(object = stem,metadata = cell.type,col.name = 'cell.type') 
stem$cc_labels<-paste0(stem@meta.data$cell.type,"-",stem@meta.data$seurat_clusters)
unique(stem$cc_labels)
#[1] "Stem-2" "Stem-0" "Stem-1"



#####################----Exporting cc_labels to tsv file


cc_labs<-stem$cc_labels
head(cc_labs)
#d21_55791_AACGGGAGTATCGAAA-1 d21_55791_AATGGAATCGCGTCGA-1 
#                    "Stem-2"                     "Stem-0" 
#d21_55791_ACGGAAGCAGTGTGCC-1 d21_55791_AGACCCGCATGGGAAC-1 
#                    "Stem-1"                     "Stem-0" 
#d21_55791_AGTACTGCAGCTGCCA-1 d21_55791_AGTCTCCGTTTATGCG-1 
#                    "Stem-2"                     "Stem-0"

cc_labs<-data.frame(cc_labs)

head(cc_labs)
#                             cc_labs
#d21_55791_AACGGGAGTATCGAAA-1  Stem-2
#d21_55791_AATGGAATCGCGTCGA-1  Stem-0
#d21_55791_ACGGAAGCAGTGTGCC-1  Stem-1
#d21_55791_AGACCCGCATGGGAAC-1  Stem-0
#d21_55791_AGTACTGCAGCTGCCA-1  Stem-2

write.table(cc_labs, file = "/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/stem.tsv", row.names=TRUE, sep="\t")

#stem.tsv
#	                        cc_labs
#d21_55791_AACGGGAGTATCGAAA-1	Stem-2
#d21_55791_AATGGAATCGCGTCGA-1	Stem-0
#d21_55791_ACGGAAGCAGTGTGCC-1	Stem-1
#d21_55791_AGACCCGCATGGGAAC-1	Stem-0



###################------Reading metadata from tsv file to assign it to seurat object



metadata<- read.table("/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/stem.tsv", header=T, sep="\t", row.names=1, check.names=F)
head(metadata)
#                             cc_labs
#d21_55791_AACGGGAGTATCGAAA-1  Stem-2
#d21_55791_AATGGAATCGCGTCGA-1  Stem-0
#d21_55791_ACGGAAGCAGTGTGCC-1  Stem-1
#d21_55791_AGACCCGCATGGGAAC-1  Stem-0
#d21_55791_AGTACTGCAGCTGCCA-1  Stem-2

stem@meta.data<- metadata
DefaultAssay(stem)<- "RNA"

Idents(stem) <- "cc_labs"

head(Idents(stem))
#d21_55791_AACGGGAGTATCGAAA-1 d21_55791_AATGGAATCGCGTCGA-1 
#                      Stem-2                       Stem-0 
#d21_55791_ACGGAAGCAGTGTGCC-1 d21_55791_AGACCCGCATGGGAAC-1 
#                      Stem-1                       Stem-0 
#d21_55791_AGTACTGCAGCTGCCA-1 d21_55791_AGTCTCCGTTTATGCG-1 
#                      Stem-2                       Stem-0 
#Levels: Stem-2 Stem-0 Stem-1


cell_labels<- c("Stem-2","Stem-0","Stem-1")



####################-------Subsetting,Normalizing and Scaling the seurat object



subsetted<- subset(stem, idents = cell_labels, invert = FALSE)
subsetted<- NormalizeData(subsetted)
subsetted<- ScaleData(subsetted)


data.input <- GetAssayData(subsetted, assay = "RNA", slot = "data")
labels <- Idents(subsetted)
head(labels)
#d21_55791_AACGGGAGTATCGAAA-1 d21_55791_AATGGAATCGCGTCGA-1 
#                      Stem-2                       Stem-0 
#d21_55791_ACGGAAGCAGTGTGCC-1 d21_55791_AGACCCGCATGGGAAC-1 
#                      Stem-1                       Stem-0 
#d21_55791_AGTACTGCAGCTGCCA-1 d21_55791_AGTCTCCGTTTATGCG-1 
#                      Stem-2                       Stem-0 
#Levels: Stem-2 Stem-0 Stem-1


meta <- data.frame(group = labels, row.names = names(labels))
head(meta)
#                              group
#d21_55791_AACGGGAGTATCGAAA-1 Stem-2
#d21_55791_AATGGAATCGCGTCGA-1 Stem-0
#d21_55791_ACGGAAGCAGTGTGCC-1 Stem-1
#d21_55791_AGACCCGCATGGGAAC-1 Stem-0
#d21_55791_AGTACTGCAGCTGCCA-1 Stem-2
#d21_55791_AGTCTCCGTTTATGCG-1 Stem-0



###################-------Cellchat object creation


cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
#The cell groups used for CellChat analysis are  Stem-2 Stem-0 Stem-1 
cellchat <- setIdent(cellchat, ident.use = "group")

groupSize <- as.numeric(table(cellchat@idents))
groupSize
#[1] 34 66 44


###################-------Setting up cell chat DB

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1887,]
which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1900,]
cellchat@DB <- CellChatDB.use
cellchat <- CellChat::subsetData(cellchat)

set.seed("1234")
plan("multisession", workers=18)
options(future.globals.maxSize = 200000*10^24)
future.seed=TRUE


###################-------Community probabbility calculation and aggregate network generation

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
table(cellchat@idents)



#Stem-2 Stem-0 Stem-1 
#    34     66     44 
#> head(cellchat@netP)
#$pathways
# [1] "LAMININ"   "THBS"      "CDH"       "APP"       "MIF"       "JAM"      
# [7] "CDH1"      "WNT"       "EGF"       "TGFb"      "GRN"       "DESMOSOME"
#[13] "OCLN"      "EPHA"      "GALECTIN"  "CADM"      "PTN"       "EPHB"     
#[19] "NOTCH"     "HSPG"      "SEMA3"     "SEMA6"     "NRG"       "PROS"     
#[25] "NCAM"      "MPZ"       "L1CAM"    

#$prob
#, , LAMININ

#          Stem-2    Stem-0      Stem-1
#Stem-2 0.2078794 0.5403887 0.009209626
#Stem-0 0.9041278 1.4314237 0.067792759
#Stem-1 0.0000000 0.0000000 0.000000000



setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/")
saveRDS(cellchat, "stem_CommunProb_before_filter.rds")

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/")
saveRDS(cellchat,"stem-aggregatenet.rds")



#############-------Cellchat object for Macropages ##############################

setwd("/home/sures105/Downloads/")
m<-readRDS("subsample_after_clustree.rds")
dim(m)
#12484   701
table(m@meta.data$seurat_clusters)

#0   1   2   3   4   5   6   7   8   9 
#155 122  86  83  64  54  54  46  25  12 


#Adding meta column:cc_labels:Eg: Macro-0, Macro-1, Macro-2....Macro-9

cell.type="Macro"
m <- AddMetaData(object = m,metadata = cell.type,col.name = 'cell.type') 
m$cc_labels<-paste0(m@meta.data$cell.type,"-",m@meta.data$seurat_clusters)
unique(m$cc_labels)
#[1] "Macro-8" "Macro-3" "Macro-9" "Macro-0" "Macro-4" "Macro-1" "Macro-6"
#[8] "Macro-7" "Macro-5" "Macro-2"


#####################----Exporting cc_labels to tsv file


cc_labs<-m$cc_labels
head(cc_labs)
#d21_55791_AAACGCTGTAAGTCAA-1 d21_55791_AAAGAACCAGCAAGAC-1 
#                   "Macro-8"                    "Macro-3" 
#d21_55791_AAAGTCCAGGGATGTC-1 d21_55791_AAGCATCCACTTCAAG-1 
#                   "Macro-9"                    "Macro-0" 

cc_labs<-data.frame(cc_labs)

head(cc_labs)
                             cc_labs
d21_55791_AAACGCTGTAAGTCAA-1 Macro-8
d21_55791_AAAGAACCAGCAAGAC-1 Macro-3
d21_55791_AAAGTCCAGGGATGTC-1 Macro-9
d21_55791_AAGCATCCACTTCAAG-1 Macro-0


write.table(cc_labs, file = "/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/macro.tsv", row.names=TRUE, sep="\t")
#macro.tsv
#	cc_labs
#d21_55791_AAACGCTGTAAGTCAA-1	Macro-8
#d21_55791_AAAGAACCAGCAAGAC-1	Macro-3
#d21_55791_AAAGTCCAGGGATGTC-1	Macro-9
#d21_55791_AAGCATCCACTTCAAG-1	Macro-0


###################------Reading metadata from tsv file to assign it to seurat object



metadata<- read.table("/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/macro.tsv", header=T, sep="\t", row.names=1, check.names=F)
head(metadata)
#                             cc_labs
#d21_55791_AAACGCTGTAAGTCAA-1 Macro-8
#d21_55791_AAAGAACCAGCAAGAC-1 Macro-3
#d21_55791_AAAGTCCAGGGATGTC-1 Macro-9
#d21_55791_AAGCATCCACTTCAAG-1 Macro-0
#d21_55791_AAGGAATGTACTGGGA-1 Macro-4
#d21_55791_AAGGTAATCCATACTT-1 Macro-1


m@meta.data<- metadata
DefaultAssay(m)<- "RNA"

Idents(m) <- "cc_labs"

unique(Idents(stem))

#[1] Macro-8 Macro-3 Macro-9 Macro-0 Macro-4 Macro-1 Macro-6 Macro-7 Macro-5 Macro-2
#10 Levels: Macro-8 Macro-3 Macro-9 Macro-0 Macro-4 Macro-1 Macro-6 ... Macro-2

cell_labels<- c("Macro-8","Macro-3","Macro-9","Macro-0","Macro-4","Macro-1","Macro-6","Macro-2","Macro-5","Macro-7")


####################-------Subsetting,Normalizing and Scaling the seurat object



subsetted<- subset(m, idents = cell_labels, invert = FALSE)
subsetted<- NormalizeData(subsetted)
subsetted<- ScaleData(subsetted)


data.input <- GetAssayData(subsetted, assay = "RNA", slot = "data")
labels <- Idents(subsetted)
head(labels)
#d21_55791_AAACGCTGTAAGTCAA-1 d21_55791_AAAGAACCAGCAAGAC-1 
#                     Macro-8                      Macro-3 
#d21_55791_AAAGTCCAGGGATGTC-1 d21_55791_AAGCATCCACTTCAAG-1 
#                     Macro-9                      Macro-0 
#d21_55791_AAGGAATGTACTGGGA-1 d21_55791_AAGGTAATCCATACTT-1 
#                     Macro-4                      Macro-1 
#10 Levels: Macro-8 Macro-3 Macro-9 Macro-0 Macro-4 Macro-1 Macro-6 ... Macro-2


meta <- data.frame(group = labels, row.names = names(labels))
head(meta)
                               group
d21_55791_AAACGCTGTAAGTCAA-1 Macro-8
d21_55791_AAAGAACCAGCAAGAC-1 Macro-3
d21_55791_AAAGTCCAGGGATGTC-1 Macro-9
d21_55791_AAGCATCCACTTCAAG-1 Macro-0


###################-------Cellchat object creation


cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
#The cell groups used for CellChat analysis are  Macro-8 Macro-3 Macro-9 Macro-0 Macro-4 Macro-1 Macro-6 Macro-7 Macro-5 Macro-2 

cellchat <- setIdent(cellchat, ident.use = "group")

groupSize <- as.numeric(table(cellchat@idents))
groupSize
#[1]  25  83  12 155  64 122  54  46  54  86


###################-------Setting up cell chat DB

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1887,]
which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1900,]
cellchat@DB <- CellChatDB.use
cellchat <- CellChat::subsetData(cellchat)

set.seed("1234")
plan("multisession", workers=18)
options(future.globals.maxSize = 200000*10^24)
future.seed=TRUE

###################-------Community probabbility calculation and aggregate network generation

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
table(cellchat@idents)

#Macro-8 Macro-3 Macro-9 Macro-0 Macro-4 Macro-1 Macro-6 Macro-7 Macro-5 Macro-2 
#     25      83      12     155      64     122      54      46      54      86 

head(cellchat@netP)
# , , HGF

#        Macro-8 Macro-3 Macro-9 Macro-0 Macro-4 Macro-1 Macro-6 Macro-7
#Macro-8       0       0       0       0       0       0       0       0
#Macro-3       0       0       0       0       0       0       0       0


setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/")
saveRDS(cellchat, "macro_CommunProb_before_filter.rds")

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
setwd("/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/")
saveRDS(cellchat,"macro-aggregatenet.rds")


#############-------Merged Cellchat object ##############################

cellchat.s <- readRDS("/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/stem-aggregatenet.rds")
cellchat.s <- updateCellChat(cellchat.s)


cellchat.m <- readRDS("/depot/tlratlif/data/LSC_scRNAseq_2023/6_cellchat/macro-aggregatenet.rds")
cellchat.m <- updateCellChat(cellchat.m)

#####Lift up CellChat objects and merge them together

group.new = levels(cellchat.s.s@idents)
cellchat.s <- liftCellChat(cellchat.s, group.new)

#The CellChat object will be lifted up using the cell labels Stem-2, Stem-0, Stem-1
#Update slots object@net, object@netP, object@idents in a single dataset...

group.new = levels(cellchat.m.m@idents)
cellchat.m <- liftCellChat(cellchat.m, group.new)

object.list <- list(Stem = cellchat.s, Macro = cellchat.m)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)


ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
