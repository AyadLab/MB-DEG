#load in libraries
library(tidyverse)
library(Seurat)
library(future)
library(patchwork)
library(ggplot2)
library(cowplot)
library(clustree)
library(ggpubr)
library(glmGamPoi)

setwd("/pathtoworkingdirectory")

options(future.globals.maxSize = +Inf)

#Load in metadata from sample series matrix file
series_matrix <- read.delim(file = "Data/GSE119926_series_matrix.txt", skip = 32)
glimpse(series_matrix)

#Optimize the titles
rownames(series_matrix) <- paste0(rownames(series_matrix), "_", series_matrix$X.Sample_title)
series_matrix$X.Sample_title <- NULL
series_mat <- as.data.frame(t(series_matrix))

#Determine the heading of series_matrix
rownames(series_mat)

#Load in Medulloblastoma Patient Tumor Samples

#file locations
file_list <- c(
  "BCH807"  = "GSM3905406_BCH807.txt",
  "BCH825"  = "GSM3905407_BCH825.txt",
  "BCH1031" = "GSM3905408_BCH1031.txt",
  "BCH1205" = "GSM3905409_BCH1205.txt",
  "MUV11"   = "GSM3905410_MUV11.txt",
  "MUV19"   = "GSM3905411_MUV19.txt",
  "MUV27"   = "GSM3905412_MUV27.txt",
  "MUV29"   = "GSM3905413_MUV29.txt",
  "MUV34"   = "GSM3905414_MUV34.txt",
  "MUV37"   = "GSM3905415_MUV37.txt",
  "MUV39"   = "GSM3905416_MUV39.txt",
  "MUV41"   = "GSM3905417_MUV41.txt",
  "MUV44"   = "GSM3905418_MUV44.txt",
  "SJ17"    = "GSM3905419_SJ17.txt",
  "SJ99"    = "GSM3905420_SJ99.txt",
  "SJ129"   = "GSM3905421_SJ129.txt",
  "SJ217"   = "GSM3905422_SJ217.txt",
  "SJ454"   = "GSM3905423_SJ454.txt",
  "SJ516"   = "GSM3905424_SJ516.txt",
  "SJ577"   = "GSM3905425_SJ577.txt",
  "SJ617"   = "GSM3905426_SJ617.txt",
  "SJ625"   = "GSM3905427_SJ625.txt",
  "SJ723"   = "GSM3905428_SJ723.txt",
  "SJ917"   = "GSM3905429_SJ917.txt",
  "SJ970"   = "GSM3905430_SJ970.txt"
)
#custom function called patient_counts creates individual seurat objects for all the patient samples and load their metadata
patient_counts <- function(sample_id, series_mat, file_list) {
  file_name <- file_list[sample_id]
  counts <- read.table(file = paste0("Data/", file_name), header = TRUE, sep = "\t", row.names = 1)
  obj <- CreateSeuratObject(counts = counts)
  obj$Methylation_Subgroup <- series_mat[sample_id, ]$`11_!Sample_characteristics_ch1`
  obj$Methylation_Subtype <- series_mat[sample_id, ]$`12_!Sample_characteristics_ch1`
  obj$Metastasis <- series_mat[sample_id, ]$`13_!Sample_characteristics_ch1`
  obj$Histology <- series_mat[sample_id, ]$`14_!Sample_characteristics_ch1`
  obj$ID <- sample_id
  return(obj)
}

#patient names are saved in the patient_id variable
patient_id <- names(file_list)

#Store all the objects in MedList
MedList <- lapply(patient_id, patient_counts, series_mat = series_mat, file_list = file_list)
names(MedList) <- patient_id
glimpse(MedList)

# SCTransform-based normalization and scaling
# SCTransform normailzation and scaling
for (i in 1:length(MedList)) {
  MedList[[i]] <- SCTransform(
    MedList[[i]],
    assay = "RNA",
    return.only.var.genes = FALSE,
    verbose = FALSE
  )
  MedList[[i]] <- RunPCA(MedList[[i]])
}

saveRDS(MedList, file = "01_split.sct_2025.02.23.RDS")

################################################################################
# MERGE AND INTEGRATION OF INDIVIDUAL OBJECTS

# merge individual patient objects
vec.objs <- c(
  MedList[[2]], MedList[[3]], MedList[[4]], MedList[[5]],
  MedList[[6]], MedList[[7]], MedList[[8]], MedList[[9]],
  MedList[[10]], MedList[[11]], MedList[[12]], MedList[[13]],
  MedList[[14]], MedList[[15]], MedList[[16]], MedList[[17]],
  MedList[[18]], MedList[[19]], MedList[[20]], MedList[[21]],
  MedList[[22]], MedList[[23]], MedList[[24]], MedList[[25]]
)
cell.IDs <- names(MedList)[1:25]
obj <- merge(
  x = MedList[[1]],
  y = vec.objs,
  add.cell.ids = cell.IDs,
  project = "Ruiz_MB",
  merge.data = TRUE
)

# # SCT layer is merged (counts, data, scale.data), RNA layer is unmerged (only counts)
# saveRDS(obj, file = "02_scaled.merge_2024.04.04.RDS")

# need to add variable features data to the merged object bc that apparently
# doesn't carry over
VariableFeatures(obj) <- c(
  VariableFeatures(MedList[[1]]),
  VariableFeatures(MedList[[2]]), VariableFeatures(MedList[[3]]), VariableFeatures(MedList[[4]]), VariableFeatures(MedList[[5]]),
  VariableFeatures(MedList[[6]]), VariableFeatures(MedList[[7]]), VariableFeatures(MedList[[8]]), VariableFeatures(MedList[[9]]),
  VariableFeatures(MedList[[10]]), VariableFeatures(MedList[[11]]), VariableFeatures(MedList[[12]]), VariableFeatures(MedList[[13]]),
  VariableFeatures(MedList[[14]]), VariableFeatures(MedList[[15]]), VariableFeatures(MedList[[16]]), VariableFeatures(MedList[[17]]),
  VariableFeatures(MedList[[18]]), VariableFeatures(MedList[[19]]), VariableFeatures(MedList[[20]]), VariableFeatures(MedList[[21]]),
  VariableFeatures(MedList[[22]]), VariableFeatures(MedList[[23]]), VariableFeatures(MedList[[24]]), VariableFeatures(MedList[[25]])
)

# while SCT assay is integrated (?), the RNA assay is not and has a counts
# layer for each patient; prep for integration

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# integrate layers of spatial assay
obj <- IntegrateLayers(
  obj,
  assay = "RNA",
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  # normalization.method = "SCT",
  verbose = FALSE, k.weight = 50 # one sample has n cells = 52
)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
# now have "integrated.cca" dim reduction in Spatial Assay

saveRDS(obj, file = "03_merged.joined_2024.04.04.RDS")

################################################################################
# PERCENT MT, CELL CYCLE SCORING

DefaultAssay(obj) <- "RNA"

rownames(obj)[grep("^MT", rownames(obj))]
# not all of these genes are mitochondrial transcripts... (e.g. MTOR...)

# mitochondrial DNA scoring
obj <- PercentageFeatureSet(
  obj,
  pattern = "^MT",
  col.name = "percent.mt",
  assay = "RNA"
)
range(obj$percent.mt)

# cell cycle scoring -- method 1 (phases independently)
## removes all signal associated with cell cycle
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

obj <- ScaleData(
  obj,
  vars.to.regress = "percent.mt"
)

obj <- CellCycleScoring(
  obj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE,
  assay = "RNA"
)

obj@meta.data$cell.cycle.phase <- Idents(obj)
Idents(obj) <- obj@meta.data$orig.ident
obj@meta.data$old.ident <- NULL
obj@meta.data$Phase <- NULL

# cell cycle scoring -- method 2 (phase difference)
## signals separating non-cycling cells and cycling cells will be maintained
## but differences in cell cycle phase among proliferating cells will be removed
obj$cell.cycle.difference <- obj$S.Score - obj$G2M.Score

saveRDS(obj, file = "04_merged.joined.scored_2025.02.23.RDS")
obj <- readRDS("04_merged.joined.scored_2025.02.23.RDS")

################################################################################
# PROCESS DARMANIS DATA

darmanis <- read.table(file = "Data/Darmanis_et_al_2015_matrix.txt", sep = "\t")
hb <- CreateSeuratObject(counts = darmanis)

DefaultAssay(hb) <- "RNA"
hb <- NormalizeData(hb)
hb <- FindVariableFeatures(hb)
hb <- ScaleData(hb)
hb <- RunPCA(hb)
hb <- FindNeighbors(hb)
hb <- FindClusters(hb, resolution = 0.5)

hb$clusters <- hb$RNA_snn_res.0.5
Idents(hb) <- hb$clusters
hb <- RenameIdents(hb, "0" = "Neural", "1" = "Neural", "2" = "Astrocytes",
                   "3" = "Endothelial", "4" = "Neural", "5" =
                     "Oligodendrocytes", "6" = "Microglia")
hb$CellTypes <- Idents(hb)
table(hb$CellTypes)

rownames(hb)[grep("^MT", rownames(hb))]
hb <- PercentageFeatureSet(
  hb,
  pattern = "^MT",
  col.name = "percent.mt",
  assay = "RNA"
)
range(hb$percent.mt)


hb <- CellCycleScoring(
  hb,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE,
  assay = "RNA"
)
hb@meta.data$cell.cycle.phase <- Idents(hb)
Idents(hb) <- hb@meta.data$orig.ident
hb@meta.data$old.ident <- NULL
hb@meta.data$Phase <- NULL
hb$cell.cycle.difference <- hb$S.Score - hb$G2M.Score

glimpse(hb@meta.data)

################################################################################
# INTEGRATE WITH DARMANIS DATA

hb <- SCTransform(
  hb,
  assay = "RNA",
  return.only.var.genes = FALSE,
  verbose = FALSE
)
hb <- RunPCA(hb)

# split obj by sample ID
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$ID)

# merge individual patient objects
cell.IDs <- c("MB", "Darmanis")
obj.merged <- merge(
  x = obj,
  y = hb,
  add.cell.ids = cell.IDs,
  project = "Ruiz_MB2",
  merge.data = TRUE
)

# need to add variable features data to the merged object bc that apparently
# doesn't carry over
VariableFeatures(obj.merged) <- c(
  VariableFeatures(obj),
  VariableFeatures(hb)
)

# while SCT assay is integrated (?), the RNA assay is not and has a counts
# layer for each patient; prep for integration

DefaultAssay(obj.merged) <- "RNA"
obj.merged <- NormalizeData(obj.merged)
obj.merged <- FindVariableFeatures(obj.merged)
obj.merged <- ScaleData(obj.merged, vars.to.regress = c("percent.mt", "G2M.Score", "S.Score"))
obj.merged <- RunPCA(obj.merged)

# integrate layers of spatial assay
obj.merged <- IntegrateLayers(
  obj.merged,
  assay = "RNA",
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  # normalization.method = "SCT",
  verbose = FALSE, k.weight = 50 # one sample has n cells = 52
)
obj.merged[["RNA"]] <- JoinLayers(obj.merged[["RNA"]])
# now have "integrated.cca" dim reduction in Spatial Assay

saveRDS(obj.merged, file = "05_darmanis.merged.joined_2025.02.23.RDS")

################################################################################
# PROCESSING AND ANNOTATIONS

table(obj.merged$CellTypes)
Idents(obj.merged) <- obj.merged$CellTypes
table(Idents(obj.merged))

obj.merged@meta.data <- obj.merged@meta.data %>%
  mutate(
    neoplasticvsnon = case_when(
      CellTypes == "Neural" ~ "Non-Neoplastic",
      CellTypes == "Astrocytes" ~ "Non-Neoplastic",
      CellTypes == "Endothelial" ~ "Non-Neoplastic",
      CellTypes == "Oligodendrocytes" ~ "Non-Neoplastic",
      CellTypes == "Microglia" ~ "Non-Neoplastic",
      is.na(CellTypes) ~ "Neoplastic"
    )
  )
table(obj.merged$neoplasticvsnon)

obj.merged@meta.data <- obj.merged@meta.data %>%
  mutate(
    cell.types = case_when(
      CellTypes == "Neural" ~ "Neural",
      CellTypes == "Astrocytes" ~ "Astrocytes",
      CellTypes == "Endothelial" ~ "Endothelial",
      CellTypes == "Oligodendrocytes" ~ "Oligodendrocytes",
      CellTypes == "Microglia" ~ "Microglia",
      is.na(CellTypes) ~ "Neoplastic"
    )
  )
table(obj.merged$cell.types)

obj.merged@meta.data <- obj.merged@meta.data %>%
  mutate(
    subgroup = case_when(
      Methylation_Subgroup == "methylation subgroup: WNT" ~ "WNT",
      Methylation_Subgroup == "methylation subgroup: Group 3" ~ "Group 3",
      Methylation_Subgroup == "methylation subgroup: Group 4" ~ "Group 4",
      Methylation_Subgroup == "methylation subgroup: SHH-adult" ~ "SHH",
      Methylation_Subgroup == "methylation subgroup: SHH-infant" ~ "SHH",
      is.na(Methylation_Subgroup) ~ "Non-Neoplastic"
    )
  )
table(obj.merged$subgroup)


obj.merged <- FindNeighbors(obj.merged, reduction = "integrated.cca", dims = 1:30)
obj.merged <- FindClusters(obj.merged, resolution = 0.5)
obj.merged <- RunUMAP(obj.merged, dims = 1:30, reduction = "integrated.cca")

obj.merged <- ScaleData(
  obj.merged, assay = "RNA",
  features = rownames(obj.merged[["RNA"]]$counts),
  vars.to.regress = c("percent.mt", "G2M.Score", "S.Score")
)

saveRDS(obj.merged, file = "06_darmanis.scaled_2025.02.23.RDS")

################################################################################
# VISUALIZATIONS

DefaultAssay(obj.merged) <- "RNA"

obj.merged <- FindVariableFeatures(obj.merged)
obj.merged <- RunPCA(obj.merged, reduction.name = "RNA.pca")
obj.merged <- FindNeighbors(obj.merged, reduction = "RNA.pca", dims = 1:30)
obj.merged <- FindClusters(obj.merged, reduction = "RNA.pca", resolution = 0.5)
obj.merged <- RunUMAP(obj.merged, reduction = "RNA.pca", dims = 1:30, reduction.name = "RNA.umap")

pdf("Plots_2025.02.23.pdf")
DimPlot(obj.merged, reduction = "umap", group.by = "neoplasticvsnon")
DimPlot(obj.merged, reduction = "umap", group.by = "cell.types")
DimPlot(obj.merged, reduction = "umap", group.by = "subgroup")
DimPlot(obj.merged, cols = c("magenta", "black"), group.by = "neoplasticvsnon")

DimPlot(obj.merged, reduction = "RNA.umap", group.by = "neoplasticvsnon")
DimPlot(obj.merged, reduction = "RNA.umap", group.by = "cell.types")
DimPlot(obj.merged, reduction = "RNA.umap", group.by = "subgroup")
DimPlot(obj.merged, cols = c("magenta", "black"), reduction = "RNA.umap", group.by = "neoplasticvsnon")

FeaturePlot(obj.merged, features = c("BAIAP2", "CDC42"), blend = T)
FeaturePlot(obj.merged, features = c("BAIAP2"), split.by = "neoplasticvsnon")
FeaturePlot(obj.merged, features = c("CDC42"), split.by = "neoplasticvsnon")
FeaturePlot(obj.merged, features = c("BAIAP2"))
FeaturePlot(obj.merged, features = c("CDC42"))
FeatureScatter(obj.merged, feature1 = "BAIAP2", feature2 = "CDC42", split.by = "neoplasticvsnon", group.by = "subgroup")

VlnPlot(obj.merged, features = "BAIAP2", group.by = "neoplasticvsnon", cols = c("magenta", "black"))
VlnPlot(obj.merged, features = "CDC42", group.by = "neoplasticvsnon", cols = c("magenta", "black"))

DotPlot(obj.merged, features = c("BAIAP2", "CDC42"), group.by = "neoplasticvsnon")
dev.off()

pdf("Plots_2_2025.02.24.pdf", height = 5, width = 10)
FeaturePlot(obj.merged, features = c("BAIAP2"), split.by = "neoplasticvsnon", reduction = "RNA.umap")
FeaturePlot(obj.merged, features = c("CDC42"), split.by = "neoplasticvsnon", reduction = "RNA.umap")
dev.off()

pdf("Plots_3_2025.02.24.pdf", height = 5, width = 5)
FeaturePlot(obj.merged, features = c("BAIAP2"), reduction = "RNA.umap")
FeaturePlot(obj.merged, features = c("CDC42"), reduction = "RNA.umap")
dev.off()

pdf("Plots_4_2025.02.24.pdf", height = 5, width = 15)
FeaturePlot(
  obj.merged,
  features = c("BAIAP2", "CDC42"),
  blend = T, cols = c("gray", "red2", "blue2"),
  reduction = "RNA.umap", order = TRUE
)
dev.off()

pdf("Plots_5_2025.02.24.pdf", height = 5, width = 7)
DimPlot(obj.merged, reduction = "RNA.umap", group.by = "ID")
dev.off()

obj.merged$cell.types <- factor(obj.merged$cell.types, levels = c("Neoplastic", "Astrocytes", "Endothelial", "Microglia", "Neural", "Oligodendrocytes"))
pdf("Plots_6_2025.02.24.pdf", height = 5, width = 5)
VlnPlot(obj.merged, features = "BAIAP2", group.by = "cell.types")
VlnPlot(obj.merged, features = "CDC42", group.by = "cell.types")
dev.off()

obj.merged$subgroup <- factor(obj.merged$subgroup, levels = c("SHH", "WNT", "Group 3", "Group 4", "Non-Neoplastic"))
pdf("Plots_7_2025.02.24.pdf", height = 5, width = 5)
VlnPlot(obj.merged, features = "BAIAP2", group.by = "subgroup")
VlnPlot(obj.merged, features = "CDC42", group.by = "subgroup")
dev.off()

pdf("Plots_8_2025.02.24.pdf", height = 5, width = 6)
DimPlot(obj.merged, reduction = "RNA.umap", group.by = "subgroup")
dev.off()

pdf("Plots_9_2025.02.24.pdf", height = 5, width = 6)
VlnPlot(obj.merged, features = "CD74", group.by = "cell.types") # microglia
VlnPlot(obj.merged, features = "GFAP", group.by = "cell.types") # astrocytes
VlnPlot(obj.merged, features = "MBP", group.by = "cell.types")  # oligodendrocytes
VlnPlot(obj.merged, features = "NEUROD6", group.by = "cell.types") # neural
VlnPlot(obj.merged, features = "ITM2A", group.by = "cell.types") # endothelial
dev.off()

# neoplastic vs. non-neoplastic differential expression
Idents(obj.merged) <- obj.merged$neoplasticvsnon
markers <- FindMarkers(obj.merged, ident.1 = "Neoplastic", ident.2 = "Non-Neoplastic", test.use = "MAST")
write.csv(markers, file = "markers.20250218.csv")

# subgroup vs. non-neoplastic differential expression
Idents(obj.merged) <- obj.merged$subgroup
SHH.markers <- FindMarkers(obj.merged, ident.1 = "SHH", ident.2 = "Non-Neoplastic", test.use = "MAST", assay = "RNA")
WNT.markers <- FindMarkers(obj.merged, ident.1 = "WNT", ident.2 = "Non-Neoplastic", test.use = "MAST", assay = "RNA")
G3.markers <- FindMarkers(obj.merged, ident.1 = "Group 3", ident.2 = "Non-Neoplastic", test.use = "MAST", assay = "RNA")
G4.markers <- FindMarkers(obj.merged, ident.1 = "Group 4", ident.2 = "Non-Neoplastic", test.use = "MAST", assay = "RNA")

write.csv(SHH.markers, file = "SHH.markers.20250219.csv")
write.csv(WNT.markers, file = "WNT.markers.20250219.csv")
write.csv(G3.markers, file = "G3.markers.20250219.csv")
write.csv(G4.markers, file = "G4.markers.20250219.csv")

deg.list <- list(SHH.markers, WNT.markers, G3.markers, G4.markers)
for (i in 1:length(deg.list)) {
  deg.list[[i]] <- deg.list[[i]][which(deg.list[[i]]$p_val_adj < 0.05 & deg.list[[i]]$avg_log2FC > 4), ]
  print(length(rownames(deg.list[[i]])))
}

# SHH.markers[which(rownames(SHH.markers) == "BAIAP2"), ]
# WNT.markers[which(rownames(WNT.markers) == "BAIAP2"), ]
# G3.markers[which(rownames(G3.markers) == "BAIAP2"), ]
# G4.markers[which(rownames(G4.markers) == "BAIAP2"), ]







