# https://colabdev.bioturing.com/notebook/identifying-tumor-cells-at-the-singlecell-level-us-cee1085d84

library(Seurat)
library(SingleR)
library(infercnv)
library(celldex)
library(ggplot2)

getwd()

# downloaded the data in : /home/tanasa/infercnv_data

# url <- 'https://cdn.bioturing.com/colab/data/singleR_infercnv_data.zip'
# download.file(url, './singleR_infercnv_data.zip')                                             # download from url
# unzip(zipfile='singleR_infercnv_data.zip', exdir='./')

# Here, we propose ikarus, a machine learning pipeline aimed at distinguishing tumor cells from normal cells at the single-cell level. 
# We test ikarus on multiple single-cell datasets, showing that it achieves high sensitivity and specificity in multiple experimental contexts.
# InferCNV is a Bayesian method, which agglomerates the expression signal of genomically adjointed genes to ascertain whether there is a gain or loss of a certain larger genomic segmen

# In this notebook, we will use a sample dataset from Wu et al. This is a single-cell atlas for human breast cancers. 
# In this notebook, in order to simplify, we will use just a sample from this atlas. 
# The chosen sample is CID4495 (GEO accession ID: GSM5354530).

# Reference: Wu, S.Z., Al-Eryani, G., Roden, D.L. et al. A single-cell and spatially resolved atlas of human breast cancers. 
# Nat Genet 53, 1334–1347 (2021). https://doi.org/10.1038/s41588-021-00911-1

sparse_matrix <- Read10X(
    data.dir='infercnv_data/GSM5354530_CID4495',
    gene.column=1,
    cell.column=1
)

seu_obj <- CreateSeuratObject(sparse_matrix)

# str(seu_obj)

# Number of cells (columns of count matrix)
n_cells <- ncol(seu_obj)

# Number of genes/features (rows of count matrix)
n_genes <- nrow(seu_obj)

# Print the result
cat("📊 Number of cells:", n_cells, "\n")
cat("🧬 Number of genes:", n_genes, "\n")

gene_names <- rownames(seu_obj)
head(gene_names)     # Show the first few
tail(gene_names)     # Show the last few
length(gene_names)   # Total number of genes

# Log normalize
seu_obj <- NormalizeData(seu_obj)
# Find highly variable genes
seu_obj <- FindVariableFeatures(seu_obj)
# Center and scale the data matrix
seu_obj <- ScaleData(seu_obj)
# Compute PCA
seu_obj <- RunPCA(seu_obj)
# Compute shared-nearest neighbor graphs
seu_obj <- FindNeighbors(seu_obj)
# Run Louvain algoritm
seu_obj <- FindClusters(seu_obj)
# Run UMAP
seu_obj <- RunUMAP(seu_obj, dims = 1:30)

head(seu_obj@meta.data,2)

options(repr.plot.width=8, repr.plot.height=8)            
DimPlot(seu_obj, group.by='seurat_clusters', label=TRUE)

print("Cell type prediction with singleR")

# 🧬 Reference Dataset Functions in SingleR

# Function	Species	Description

# BlueprintEncodeData()	Human	Combines Blueprint and ENCODE; immune cells and progenitors from various tissues. Good general-purpose immune+hematopoietic reference.
# DatabaseImmuneCellExpressionData()	Human	Aggregated immune reference built from multiple datasets (e.g., Monaco, Blueprint). Useful for robust immune cell annotation.
# HumanPrimaryCellAtlasData()	Human	Diverse primary human cells across tissues. Great for broad tissue-level annotation.
# ImmGenData()	Mouse	ImmGen consortium dataset; high-quality mouse immune cell reference.
# MonacoImmuneData()	Human	Immune cells sorted from peripheral blood (PBMCs); good for blood-based studies.
# MouseRNAseqData()	Mouse	General mouse tissue/cell atlas (Tabula Muris).
# NovershternHematopoieticData()	Human	Bone marrow and hematopoietic stem/progenitor cells. Deep for myeloid/erythroid lineages.

# Cell type prediction with singleR 

# We will apply automatical annotation, which is a method that can predict cell type automatically without the need of user’s experience. 
# singleR is an automatically labelling method that utilizes annotated reference data of pur cell types to infer the cluster labels 
# in a new dataset. 

# The available functions as well as the available datasets are:

# BlueprintEncodeData
# DatabaseImmuneCellExpressionData
# HumanPrimaryCellAtlasData
# ImmGenData
# MonacoImmuneData
# MouseRNAseqData
# NovershternHematopoieticData

# When calling these functions, an important parameter is:
# ensembl: boolean, TRUE if using Ensembl gene ID, FALSE if using Gene Symbol

# In this notebook, we will use the Human Primary Cell Atlas Data as reference. 
# Our dataset genes are gene symbols so we will set ensembl=FALSE

# celldex

hpca <- celldex::HumanPrimaryCellAtlasData(ensembl=FALSE)

hpca
head(rownames(hpca),2)
head(colnames(hpca),2)
colData(hpca)
rowData(hpca)

print("The cell types in hpca:")

unique(colData(hpca)$label.fine)
length(unique(colData(hpca)$label.fine))

print("The cell types in hpca: the main types")

unique(colData(hpca)$label.main)
length(unique(colData(hpca)$label.main))



sce_obj <- as.SingleCellExperiment(seu_obj)

sce_obj
head(rownames(sce_obj),2)
head(colnames(sce_obj),2)
colData(sce_obj)
rowData(sce_obj)

# str(sce_obj)

print("There are two options to run singleR")

print("1. predict cell type for individual cells")

print("2. predict cell type for clusters")



print("Run singleR for individual cells")

# In order to run singleR for individual cells, simply call SingleR and input the 3 mandatory parameters:

# test: our input dataset that we want to annotate
# ref: reference imported from celldex
# labels: reference label, stored in label.main slot in the reference data. E.g.: REFERENCE_DATA$label.main

singler_result <- SingleR(
    test = sce_obj,
    ref = hpca,
    labels = hpca$label.main,
)

unique(singler_result$pruned.labels)
unique(singler_result$labels)

head(data.frame(singler_result))
table(singler_result$pruned.labels)

# add to metadata and visualize
seu_obj@meta.data$singler_cell <- singler_result$pruned.labels

options(repr.plot.width=12, repr.plot.height=8)             # figure size
DimPlot(seu_obj, group.by='singler_cell')

print("Run singleR on clusters of cells")

# A second option to run singleR is predicting cell types for the whole cluster. 
# This option required the predefined cluster labels. 

louvain_clusters <- seu_obj@meta.data$seurat_clusters

# In order to run singleR for cluster, we need to add the parameter clusters. 
# This parameter will denote the predefined clusters.

singler_result2 <- SingleR(
    test = sce_obj,
    ref = hpca,
    labels = hpca$label.main,
    clusters = louvain_clusters,
)

singler_result2 

unique(singler_result2$labels)
unique(singler_result2$pruned.labels)

# The singleR result will include prediction scores of all available cell types for each clusters. 
# The final prediction will be in the pruned.labels column.

head(data.frame(singler_result2))

singler_result2$pruned.labels

# Some clusters are annotated as the same cell type, we will make these annotations unique to keep the predefined clusters.

singler_result_unique <- list()
for (i in seq(1:length(singler_result2$pruned.labels))) {
    singler_result_unique <- append(
        singler_result_unique, 
        paste(as.character(i-1), singler_result2$pruned.labels[i])
    )
}
singler_result_unique <- unlist(singler_result_unique)

singler_result_unique

# add to meta data
seu_obj@meta.data$singler_cluster_unique <- seu_obj@meta.data$seurat_cluster
levels(seu_obj@meta.data$singler_cluster_unique) <- singler_result_unique

# Set figure size for side-by-side layout
options(repr.plot.width = 18, repr.plot.height = 8)

# DimPlot colored by unique SingleR cluster annotation
p1 <- DimPlot(seu_obj, group.by = "singler_cluster_unique", label = TRUE) +
  ggtitle("SingleR: Unique Cluster Annotation")

# DimPlot colored by original Seurat clustering
p2 <- DimPlot(seu_obj, group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Seurat Clusters")

# Combine and display side by side
p1 + p2




# Some clusters are pretty similar such as 0 T_cells, 1 T_cells, these clusters can be merged together.

# However, 6 T cells and 7 T_cells clusters are extremely different from other T cells clusters, 
# and they seem to be clustered with endothelial cells. 
# Note that, breast cancer tumor cell are also originated from epithelial cells. 
# So these clusters of T cells are likely misclassified (because there are no tumor cell type in the dataset). 
# We will confirm our hypothesis with inferCNV in the following.
# Now, we will merge the similar clusters.

singler_cluster <- c(
    'T_cells_1',
    'T_cells_1',
    'B_cells_1',
    'B_cells_2',
    'B_cells_3',
    'Macrophage',
    'T_cells_2',
    'T_cells_2',
    'Monocytes',
    'T_cells_3',
    'Epithelial_cells',
    'Fibroblasts',
    'Tissue_stem_cells',
    'T_cells_3',
    'Endothelial_cells',
    'Epithelial_cells',
    'Monocytes'
)

seu_obj@meta.data$singler_cluster <- seu_obj@meta.data$seurat_cluster
levels(seu_obj@meta.data$singler_cluster) <- singler_cluster

unique(seu_obj@meta.data$seurat_cluster)
unique(seu_obj@meta.data$singler_cluster)

library("ggplot2")
options(repr.plot.width = 18, repr.plot.height = 8) 

p1 <- DimPlot(seu_obj, group.by = "singler_cluster", label = TRUE) +
  ggtitle("Clusters by SingleR")

p2 <- DimPlot(seu_obj, group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters by Seurat")

p1 + p2  # patchwork handles side-by-side layout


getwd()
# list.files()

singler_labels <- dplyr::select(seu_obj@meta.data, singler_cluster)
write.table(singler_labels,'./infercnv_data/singler_labels.txt' ,sep='\t', quote=FALSE, col.names=FALSE)

head(read.table('./infercnv_data/singler_labels.txt'), 3)

head(read.table('infercnv_data/gencode_v19_gene_pos.txt'), 3)

reference_cells <- c(
    'T_cells_1',
    'T_cells_3',
    'B_cells_1',
    'B_cells_2',
    'B_cells_3',
    'Macrophage',
    'Monocytes',
    'Fibroblasts',
    'Tissue_stem_cells',
    'Endothelial_cells'
)

# prepare the input for inferCNV

# Get counts matrix directly
counts_matrix <- GetAssayData(seu_obj, assay = "RNA", slot = "counts")

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = "/home/tanasa/infercnv_data/singler_labels.txt",
  delim = "\t",
  gene_order_file = "/home/tanasa/infercnv_data/gencode_v19_gene_pos.txt",
  ref_group_names = reference_cells
)

str(infercnv_obj) 

options(scipen=100)
infercnv_obj = infercnv::run(
    infercnv_obj = infercnv_obj,
    cutoff = 0.1,                            # we are using 10x data
    out_dir = "infercnv_data",  
    cluster_by_groups = TRUE, 
    denoise = TRUE,
    HMM = TRUE
)

str(infercnv_obj)

getwd()


