#load libraries
library(Seurat)
library(tidyverse)
library(sctransform)
library(BiocManager)
library(SoupX)

#QC Measures -- Removing ambient RNA contamination with SoupX - https://doi.org/10.1093/gigascience/giaa151
#Estimated global rho values: P6 = 0.03; P7 = 0.01; P13 = 0.03; P15 = 0.01; P16 = 0.03; P20 = 0.02

P6_sc <- load10X('~/placenta_cellranger_outs/P6/')
P6 <- Read10X(data.dir = "~/placenta_cellranger_outs/P6/filtered_feature_bc_matrix/")
P6 <- CreateSeuratObject(counts = P6)
P6  <- SCTransform(P6, verbose = F)
P6  <- RunPCA(P6, verbose = F)
P6  <- RunUMAP(P6, dims = 1:20, verbose = F)
P6  <- FindNeighbors(P6, dims = 1:20, verbose = F)
P6  <- FindClusters(P6, verbose = T)
P6_meta <- P6@meta.data
P6_sc <- setClusters(P6_sc, setNames(P6_meta$seurat_clusters, rownames(P6_meta)))
P6_sc <- autoEstCont(P6_sc)
P6_out <- adjustCounts(P6_sc)

P7_sc <- load10X('~/placenta_cellranger_outs/P7/')
P7 <- Read10X(data.dir = "~/placenta_cellranger_outs/P7/filtered_feature_bc_matrix/")
P7 <- CreateSeuratObject(counts = P7)
P7 <- SCTransform(P7, verbose = F)
P7 <- RunPCA(P7, verbose = F)
P7 <- RunUMAP(P7, dims = 1:20, verbose = F)
P7 <- FindNeighbors(P7, dims = 1:20, verbose = F)
P7 <- FindClusters(P7, verbose = T)
P7_meta <- P7@meta.data
P7_sc <- setClusters(P7_sc, setNames(P7_meta$seurat_clusters, rownames(P7_meta)))
P7_sc <- autoEstCont(P7_sc)
P7_out <- adjustCounts(P7_sc)

P13_sc <- load10X('~/placenta_cellranger_outs/P13/')
P13 <- Read10X(data.dir = "~/placenta_cellranger_outs/P13/filtered_feature_bc_matrix/")
P13 <- CreateSeuratObject(counts = P13)
P13 <- SCTransform(P13, verbose = F)
P13 <- RunPCA(P13, verbose = F)
P13 <- RunUMAP(P13, dims = 1:20, verbose = F)
P13 <- FindNeighbors(P13, dims = 1:20, verbose = F)
P13 <- FindClusters(P13, verbose = T)
P13_meta <- P13@meta.data
P13_sc <- setClusters(P13_sc, setNames(P13_meta$seurat_clusters, rownames(P13_meta)))
P13_sc <- autoEstCont(P13_sc)
P13_out <- adjustCounts(P13_sc)

P15_sc <- load10X('~/placenta_cellranger_outs/P15/')
P15 <- Read10X(data.dir = "~/placenta_cellranger_outs/P15/filtered_feature_bc_matrix/")
P15 <- CreateSeuratObject(counts = P15)
P15 <- SCTransform(P15, verbose = F)
P15 <- RunPCA(P15, verbose = F)
P15 <- RunUMAP(P15, dims = 1:20, verbose = F)
P15 <- FindNeighbors(P15, dims = 1:20, verbose = F)
P15 <- FindClusters(P15, verbose = T)
P15_meta <- P15@meta.data
P15_sc <- setClusters(P15_sc, setNames(P15_meta$seurat_clusters, rownames(P15_meta)))
P15_sc <- autoEstCont(P15_sc)
P15_out <- adjustCounts(P15_sc)

P16_sc <- load10X('~/placenta_cellranger_outs/P16/')
P16 <- Read10X(data.dir = "~/placenta_cellranger_outs/P16/filtered_feature_bc_matrix/")
P16 <- CreateSeuratObject(counts = P16)
P16 <- SCTransform(P16, verbose = F)
P16 <- RunPCA(P16, verbose = F)
P16 <- RunUMAP(P16, dims = 1:20, verbose = F)
P16 <- FindNeighbors(P16, dims = 1:20, verbose = F)
P16 <- FindClusters(P16, verbose = T)
P16_meta <- P16@meta.data
P16_sc <- setClusters(P16_sc, setNames(P16_meta$seurat_clusters, rownames(P16_meta)))
P16_sc <- autoEstCont(P16_sc)
P16_out <- adjustCounts(P16_sc)

P20_sc <- load10X('~/placenta_cellranger_outs/P20/')
P20 <- Read10X(data.dir = "~/placenta_cellranger_outs/P20/filtered_feature_bc_matrix/")
P20 <- CreateSeuratObject(counts = P20)
P20 <- SCTransform(P20, verbose = F)
P20 <- RunPCA(P20, verbose = F)
P20 <- RunUMAP(P20, dims = 1:20, verbose = F)
P20 <- FindNeighbors(P20, dims = 1:20, verbose = F)
P20 <- FindClusters(P20, verbose = T)
P20_meta <- P20@meta.data
P20_sc <- setClusters(P20_sc, setNames(P20_meta$seurat_clusters, rownames(P20_meta)))
P20_sc <- autoEstCont(P20_sc)
P20_out <- adjustCounts(P20_sc)

#Creating Seurat Objects with the SoupX filtered data
P6_Control_One <- CreateSeuratObject(counts = P6_out, min.cells = 5, min.features = 200)
P7_Stress_One <- CreateSeuratObject(counts = P7_out, min.cells = 5, min.features = 200)
P13_Stress_Two <- CreateSeuratObject(counts = P13_out, min.cells = 5, min.features = 200)
P15_Stress_Two <- CreateSeuratObject(counts = P15_out, min.cells = 5, min.features = 200)
P16_Control_Two <- CreateSeuratObject(counts = P16_out, min.cells = 5, min.features = 200)
P20_Control_Two <- CreateSeuratObject(counts = P20_out, min.cells = 5, min.features = 200)

#To simplify quality control we can merge the data sets (This is NOT the same as integrating the data)
#While you can do this for each individual sample in this data, merging simplifies the process with larger data sets
rm(P6, P7, P13, P15, P16, P20, P6_out, P7_out, P13_out, P15_out, P16_out, P20_out,
   P6_sc, P7_sc, P13_sc, P15_sc, P16_sc, P20_sc, P6_meta, P7_meta, P13_meta, P15_meta,
   P16_meta, P20_meta)
merged_seurat <- merge(P13_Stress_Two, y= c(P15_Stress_Two, P16_Control_Two,
                                            P20_Control_Two, P6_Control_One,
                                            P7_Stress_One), 
                       add.cell.ids = ls()[1:6])

#Create a 'Sample' Column in the meta data (you can view the metadata first)
view(merged_seurat@meta.data)

merged_seurat$sample <- rownames(merged_seurat@meta.data)

#Separate the sample column into the sample ID, treatment type, and barcode
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample',
                                    into = c('Sample', 'Treatment','Batch', 'Barcode'),
                                    sep = '_')

#To make sure you have data from all of the samples
unique(merged_seurat@meta.data$Sample)
unique(merged_seurat@meta.data$Treatment)

#Calculate mitochondrial percentage and add it to the metadata in a new column
merged_seurat$MitoPercent <- PercentageFeatureSet(merged_seurat, pattern = '^mt-')
merged_seurat$HbPercent <- PercentageFeatureSet(merged_seurat, pattern = '^Hb')
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

#Filter out unwanted data
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 200 &
                                   nCount_RNA > 500 & nFeature_RNA < 10000 &
                                   nCount_RNA < 150000 & MitoPercent < 10 & HbPercent < 5)

#You can see how many cells were filtered out: 54392 cells to 23353 cells
merged_seurat
merged_seurat_filtered

#Visualizing the Data

FeatureScatter(merged_seurat, feature1 = 'HbPercent', feature2 = 'nCount_RNA')
FeatureScatter(merged_seurat, feature1 = 'HbPercent', feature2 = 'nFeature_RNA')
FeatureScatter(merged_seurat, feature1 = 'HbPercent', feature2 = 'log10GenesPerUMI')

metadata <- merged_seurat_filtered@meta.data

CountDensity_NoRBC <- metadata %>% 
  ggplot(aes(x=nCount_RNA)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
#visualize unique features/genes per cell
featurespercell_NoRBC <- metadata %>% 
  ggplot(aes(x=nFeature_RNA)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 250)
#mitochondrial information with features and UMIs...can clearly see RBC at bottom in a cluster (pulling the entire regression down) and can also see the cells with high mito fractions
Mito_with_features_UMI_NoRBC <- metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MitoPercent)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)
#mito ratio
metadata %>% 
  ggplot(aes(x=MitoPercent)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
#cellular complexity
Complexity_NoRBC <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

#Separate out the samples into individual Seurat objects again
# Split the Seurat object by the 'Sample' metadata column and assign each sample to a new object
split.obj <- SplitObject(merged_seurat_filtered, split.by = "Sample")

P6 <- split.obj$P6
P7 <- split.obj$P7
P13 <- split.obj$P13
P15 <- split.obj$P15
P16 <- split.obj$P16
P20 <- split.obj$P20

#Removing objects no longer needed to declutter
rm(P13_Stress_Two, P15_Stress_Two, P16_Control_Two, P20_Control_Two, 
   P6_Control_One, P7_Stress_One,
   split.obj, merged_seurat, merged_seurat_filtered)



#Transforming Data for each sample using sctransform (https://doi.org/10.1186/s13059-019-1874-1)

BiocManager::install("glmGamPoi")

P6 <- SCTransform(P6, method = "glmGamPoi", vars.to.regress = "MitoPercent", verbose = FALSE)
P7 <- SCTransform(P7, method = "glmGamPoi", vars.to.regress = "MitoPercent", verbose = FALSE)
P13 <- SCTransform(P13, method = "glmGamPoi", vars.to.regress = "MitoPercent", verbose = FALSE)
P15 <- SCTransform(P15, method = "glmGamPoi", vars.to.regress = "MitoPercent", verbose = FALSE)
P16 <- SCTransform(P16, method = "glmGamPoi", vars.to.regress = "MitoPercent", verbose = FALSE)
P20 <- SCTransform(P20, method = "glmGamPoi", vars.to.regress = "MitoPercent", verbose = FALSE)

#Perform Integration
list <- c(P6, P7, P13, P15, P16, P20)
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features)
placenta <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
saveRDS(placenta, file = "./integrated_placenta(12_16_22).rds")

#Perform dimension reduction on Integrated Data Set
#Make sure the DefaultAssay is 'integrated'

placenta <- RunPCA(object = placenta)
ElbowPlot(placenta)
placenta <- RunUMAP(object = placenta, reduction = "pca", dims = 1:20)
placenta <- FindNeighbors(placenta, reduction = "pca", dims = 1:20)
#I ran the FindCluster function at resolutions 0.1, 0.2, 0.3, 0.4, 0.5, & 0.7 and found 0.4 to best represent the data
placenta <- FindClusters(placenta, resolution = 0.4)
placenta <- RunTSNE(placenta, dims = 1:20)

#view the clusters and also look at them by treatment
clusters <- DimPlot(placenta, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
treatment <- DimPlot(placenta, reduction = 'umap', group.by = 'Treatment')
treatment|clusters

#Finding cluster IDs
DefaultAssay(placenta) <- "RNA"

# Normalize RNA data for visualization purposes
placenta_RNA_Assay <- NormalizeData(placenta, verbose = FALSE)

# T Cells (Cd3d, Cd8a, Cd4) = Cluster 12
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Cd3d', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Cd8a', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Cd4', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

# Monocytes (Ccr2)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Ly6c1', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Cx3cr1', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Cd43', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

# B cells (Cd19, Ms4a1; aka Cd20, Cd22) = Cluster 1
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Cd19', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Ms4a1', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Cd22', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)


FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Itgam', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

#Endothelial Cells (Lyve1, Tek, Kdr; aka Vegfr2) = Cluster 6
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Lyve1', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Tek', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Kdr', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)



FeaturePlot(placenta_RNA_Assay, reduction = "umap", features = 'Psg18', 
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
