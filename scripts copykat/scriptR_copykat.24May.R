
# copy from GCP : scp -i ~/.ssh/google_compute_engine -r tanasa@34.72.93.86:/home/tanasa/CNS/rGBM/seu*plots .
# rsync -ruvh -e "ssh -i ~/.ssh/google_compute_engine" tanasa@34.72.93.86:/home/tanasa/CNS/rGBM/seu*plots .

#####################################################################
#####################################################################
# Run CopyKAT CNV prediction (on single samples, RPCA / CCA / Harmony integrated rds objects)
#####################################################################
#####################################################################

library(copykat)
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(gplots)
library(RColorBrewer)
library(parallelDist)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(gplots)        # for heatmap.3
library(parallelDist)  # fast distance function
library(ComplexHeatmap)
library(proxy)  # needed for cosine
library(RColorBrewer)
library(grid)
library(svpluscnv)
# library(Cairo)  

#####################################################################
#####################################################################

library(BiocParallel)
param <- MulticoreParam(workers = 10)
register(param)

library(future)
options(future.globals.maxSize = 1000e9)
 
#####################################################################
#####################################################################

obj = readRDS("seurat5_integrated_RPCA.rds")
# head(obj@meta.data)
# tail(obj@meta.data)
# str(obj@meta.data) 

#####################################################################
#####################################################################

illness = "LGG"

# Create output directory
output_dir <- paste0(illness, ".copykat.output_22may_rds", sep="")
dir.create(output_dir, showWarnings = FALSE)

#####################################################################
#####################################################################
# 1. Extract raw counts matrix

# counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts") %>% as.matrix() # not an option in Seurat5

# This retrieves the full raw count matrix from the "counts" slot of the "RNA" assay for the entire Seurat object obj.
# Key points:
# It does not return only a single sample or layer.
# It gives you the aggregated gene-by-cell matrix (raw UMI counts) for all cells currently stored in the object, regardless of how many samples or batches were integrated into it.
# If obj contains multiple samples (e.g. from integration or merge), the returned matrix will include all cells from all samples, typically with cell names labeled accordingly (e.g. "Sample1_Cell1", "Sample2_Cell99", etc.).

#####################################################################
#####################################################################

obj@active.assay   # e.g. "RNA"
DefaultAssay(obj)  # returns "RNA"

names(obj@commands)
Reductions(obj) 

table(obj@meta.data$seurat_clusters)
table(obj@meta.data$rpca_clusters)

# table(obj@meta.data$cca_clusters)
# table(obj@meta.data$harmony_clusters)
# levels(obj@active.ident)
# str(obj)

#####################################################################
##################################################################### Join Layers

obj_copykat = obj
obj_copykat <- JoinLayers(obj_copykat)

# str(obj_copykat)
Layers(obj_copykat)
Reductions(obj_copykat) 

# names(obj[["pca"]])
# or access directly
# obj[["pca"]]@cell.embeddings           # matrix of PC coordinates per cell
# obj[["pca"]]@feature.loadings          # gene loadings per PC
# obj[["pca"]]@stdev                     # standard deviation per PC
# obj[["pca"]]@key                       # e.g. "PC_"

# Layers(obj[["RNA"]])
# [1] "counts.LGG-04-1"     "counts.LGG-04-2"     "counts.LGG-04-3"    
# [4] "counts.LGG-03"       "data.LGG-04-1"       "scale.data.LGG-04-1"
# [7] "data.LGG-04-2"       "scale.data.LGG-04-2" "data.LGG-04-3"      
#[10] "scale.data.LGG-04-3" "data.LGG-03"         "scale.data.LGG-03"  
#[13] "scale.data"  

# Layers(obj_copykat[["RNA"]])
#[1] "data"                "counts"              "scale.data.LGG-04-1"
#[4] "scale.data.LGG-04-2" "scale.data.LGG-04-3" "scale.data.LGG-03"  
#[7] "scale.data"  

# For the original object, compute the size of matrices
cat("Dimensions of layers in obj:\n")
for (layer_name in Layers(obj[["RNA"]])) {
  dims <- dim(obj[["RNA"]]@layers[[layer_name]])
  cat(sprintf("  %-25s : %d x %d\n", layer_name, dims[1], dims[2]))
}

# For the copykat-joined object
cat("\nDimensions of layers in obj_copykat:\n")
for (layer_name in Layers(obj_copykat[["RNA"]])) {
  dims <- dim(obj_copykat[["RNA"]]@layers[[layer_name]])
  cat(sprintf("  %-25s : %d x %d\n", layer_name, dims[1], dims[2]))
}

#####################################################################
#####################################################################

# head(obj_copykat@meta.data)
#                                     orig.ident nCount_RNA nFeature_RNA
#LGG-04-1_LGG-04-1_AAACCCAAGTCCCAGC-1   LGG-04-1      10672         2982
#LGG-04-1_LGG-04-1_AAACCCAAGTGCCGAA-1   LGG-04-1       1781          874
#LGG-04-1_LGG-04-1_AAACCCACAAAGTGTA-1   LGG-04-1      13339         4389
#LGG-04-1_LGG-04-1_AAACCCACACTGGAAG-1   LGG-04-1       7081         2933
#LGG-04-1_LGG-04-1_AAACCCATCGACCATA-1   LGG-04-1       2554         1596
#LGG-04-1_LGG-04-1_AAACCCATCTTCCCAG-1   LGG-04-1       5215         1859
#                                       sample percent.mt RNA_snn_res.0.5
#LGG-04-1_LGG-04-1_AAACCCAAGTCCCAGC-1 LGG-04-1  8.7706147               0
#LGG-04-1_LGG-04-1_AAACCCAAGTGCCGAA-1 LGG-04-1  0.9545199               0
#LGG-04-1_LGG-04-1_AAACCCACAAAGTGTA-1 LGG-04-1  8.6138391               3
#LGG-04-1_LGG-04-1_AAACCCACACTGGAAG-1 LGG-04-1  9.2783505              10
#LGG-04-1_LGG-04-1_AAACCCATCGACCATA-1 LGG-04-1  2.0751762               7
#LGG-04-1_LGG-04-1_AAACCCATCTTCCCAG-1 LGG-04-1  7.8044104               0
#                                     seurat_clusters rpca_clusters
#LGG-04-1_LGG-04-1_AAACCCAAGTCCCAGC-1               0             0
#LGG-04-1_LGG-04-1_AAACCCAAGTGCCGAA-1               8             8
#LGG-04-1_LGG-04-1_AAACCCACAAAGTGTA-1               2             2
#LGG-04-1_LGG-04-1_AAACCCACACTGGAAG-1               8             8
#LGG-04-1_LGG-04-1_AAACCCATCGACCATA-1               6             6
#LGG-04-1_LGG-04-1_AAACCCATCTTCCCAG-1               1             1

# dim(GetAssayData(obj_copykat, slot="counts"))
# [1] 36601 11055
# dim(GetAssayData(obj_copykat, slot="data"))
# [1] 36601 11055

dim(LayerData(obj_copykat, layer="counts"))
dim(LayerData(obj_copykat, layer="data"))

# Save raw counts matrix
write.table(
  as.matrix(LayerData(obj_copykat, layer = "counts")),
  file = file.path(output_dir, "counts_matrix.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Save normalized matrix
write.table(
  as.matrix(LayerData(obj_copykat, layer = "data")),
  file = file.path(output_dir, "norm_matrix.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

#####################################################################
#####################################################################
#####################################################################
################################################# SUBSAMPLE 1000 cells

# set.seed(42)  # For reproducibility

# harmony_mat <- obj_copykat[["harmony"]]@cell.embeddings
# head(harmony_mat)               # Rows = cells, Columns = harmony_1, harmony_2, ...
# write.table(harmony_mat, file = "harmony_embeddings.tsv", sep = "\t", quote = FALSE)

# n_total <- 1000                          # Total number of cells to sample
# clusters <- obj_copykat$harmony_clusters # Get cluster identities

# Compute number of cells per cluster (proportional)
# cluster_counts <- table(clusters)
# n_per_cluster <- round(n_total * cluster_counts / sum(cluster_counts))

# Fix rounding mismatch (sum may not equal 1000)
# while (sum(n_per_cluster) != n_total) {
#  diff <- n_total - sum(n_per_cluster)
#  max_cluster <- names(which.max(cluster_counts))
#  n_per_cluster[max_cluster] <- n_per_cluster[max_cluster] + diff
# }

# Sample cells
# cells_to_keep <- unlist(lapply(names(n_per_cluster), function(clust) {
#  cells <- WhichCells(obj_copykat, idents = clust)
#  sample(cells, min(n_per_cluster[clust], length(cells)))
# }))

# Subset the object
# obj_sub <- subset(obj_copykat, cells = cells_to_keep)
# head(obj_sub@meta.data)
# dim(obj_sub@meta.data)

# Save the object as an RDS file in output_dir
# saveRDS(obj_sub, file = file.path(output_dir, "copykat_obj_sub.rds"))

#####################################################################
##################################################################### restart

# readRDS("copykat_obj_sub.rds")
# obj_sub = readRDS("copykat_obj_sub.rds")

#####################################################################
#####################################################################

# 2. Run CopyKAT without known normal cells

#     copykat(
#       rawmat = rawdata,
#       id.type = "S",
#       cell.line = "no",
#       ngene.chr = 5,
#       LOW.DR = 0.05,
#       UP.DR = 0.1,
#       win.size = 25,
#       norm.cell.names = "",
#       KS.cut = 0.1,
#       sam.name = "",
#       distance = "euclidean",
#       output.seg = "FALSE",
#       plot.genes = "TRUE",
#       genome = "hg38",
#       n.cores = 20
#     )

# rawmat: raw data matrix; genes in rows; cell names in columns.
# id.type: gene id type: Symbol or Ensemble.
# cell.line: if the data are from pure cell line,put "yes"; if cell line
#            data are a mixture of tumor and normal cells, still put "no".
# ngene.chr: minimal number of genes per chromosome for cell filtering.
# LOW.DR: minimal population fractions of genes for smoothing.
# UP.DR: minimal population fractions of genes for segmentation.
# win.size: minimal window sizes for segmentation.
# norm.cell.names: a vector of normal cell names.
# KS.cut: segmentation parameters, input 0 to 1; larger looser criteria.
# sam.name: sample name.
# distance: distance methods include euclidean, and correlation converted
#           distance include pearson and spearman.
# output.seg: TRUE or FALSE, output seg file for IGV visualization
# plot.genes: TRUE or FALSE, output heatmap of CNVs with genename labels
# genome: hg38 or mm10, current version only work for human or mouse genes
# n.cores: number of cores for parallel computing.

# obj_copykat = obj_sub # to reuse the code below.

dim(obj_copykat@meta.data)

cat("Running CopyKAT on integrated object...\n")

# Extract merged counts matrix for CopyKAT
counts_matrix <- as.matrix(LayerData(obj_copykat[["RNA"]], layer = "counts"))

dim(counts_matrix)

# Run CopyKAT
copykat_result <- copykat(
  rawmat          = counts_matrix,
  id.type         = "S",            # "S" = gene symbols; change to "E" if using Ensembl
  cell.line = "no",
  ngene.chr       = 5,
  LOW.DR = 0.05,
  UP.DR = 0.1,
  win.size        = 25,
  norm.cell.names = NULL,           # No normal reference cells
  KS.cut = 0.1,
  distance        = "euclidean",
  output.seg = TRUE,
  plot.genes = FALSE,
  genome = "hg20",                  # https://github.com/navinlabcode/copykat/issues/4, or "hg38", "GRCh38" (all acceptable)
  n.cores         = 10,
  sam.name = file.path(output_dir, "LGG")                   # the sample name, and the output folder
  # smooth.output = TRUE,
  # sam.name        = ""
)

# Hg20, hg38, or GRCh38 are same, referring to the current reference genome version. 
# https://github.com/navinlabcode/copykat/tree/master/man

# [1] "running copykat v1.1.0"
# [1] "step1: read and filter data ..."
# [1] "36601 genes, 1000 cells in raw data"
# [1] "10821 genes past LOW.DR filtering"
# [1] "step 2: annotations gene coordinates ..."
# [1] "start annotation ..."
# [1] "step 3: smoothing data with dlm ..."
# [1] "step 4: measuring baselines ..."
# [1] "step 5: segmentation..."
# [1] "step 6: convert to genomic bins..."
# [1] "step 7: adjust baseline ..."
# [1] "step 8: final prediction ..."
# [1] "step 9: saving results..."
# [1] "step 10: ploting heatmap ...

#####################################################################
#####################################################################
# Step 3: Save CopyKAT results
#####################################################################
#####################################################################

# Write predictions table to output_dir

write.table(
  copykat_result$prediction,
  file = file.path(output_dir, "copykat_predictions.tsv"),   # copykat_harmony.tsv
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Write matrix table to output_dir

write.table(
  copykat_result$CNAmat,
  file = file.path(output_dir, "copykat_cna_matrix.tsv"),   # copykat_harmony.tsv
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Output saved to copykat_output/.\n")

names(copykat_result)
dim(copykat_result$prediction)
dim(copykat_result$CNAmat)

head(copykat_result$CNAmat, 2)
# LGG.03_LGG.03_AGGTCTAGTTCTTGCC.1 LGG.03_LGG.03_AGTACCAAGATGCCGA.1
# 1                     0.0003138822                        0.0537936
# 2                     0.0003138822                        0.0537936

saveRDS(copykat_result, file = file.path(output_dir, "copykat_runs_results.rds"))

#####################################################################
##################################################################### add new column to meta.data

# Extract only the copykat.pred column as a named vector
copykat_vector <- copykat_result$prediction$copykat.pred
names(copykat_vector) <- copykat_result$prediction$cell.names

# Add this vector to Seurat metadata
obj_copykat <- AddMetaData(obj_copykat, metadata = copykat_vector, col.name = "copykat.pred")

head(obj_copykat@meta.data$copykat.pred, 2)
table(obj_copykat@meta.data$copykat.pred)

# A note :
# obj_copykat <- AddMetaData(obj_copykat, metadata = copykat_labels, col.name = "copykat.pred")
# This works only if:
# copykat_labels is a named vector — where names are cell barcodes, 
# and values are the corresponding metadata entries.

# But in your case, copykat_labels is a data frame, not a vector. So what happens?
# AddMetaData() gets confused unless it’s a single-column named vector 
# (or a properly indexed data frame with matching rownames).

# Compare rownames :
# all(rownames(copykat_result$prediction) == rownames(obj_copykat@meta.data))
# To see if they contain the same rownames but possibly in different order:
# setequal(rownames(copykat_result$prediction), rownames(obj_copykat@meta.data))
# To find which rownames are missing or different, you can do:
# setdiff(rownames(copykat_result$prediction), rownames(obj_copykat@meta.data))
# setdiff(rownames(obj_copykat@meta.data), rownames(copykat_result$prediction))

# View(copykat_result$prediction)       # open prediction table
# View(copykat_result$CNAmat)           # view CNA matrix
# str(copykat_result$hclustering, 1)    # quick look at hclust info

# head(obj_copykat@meta.data, 2)
#                                     orig.ident nCount_RNA nFeature_RNA
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1   LGG-04-1        517          265
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1   LGG-04-1      16196         3761
#                                        sample percent.mt RNA_snn_res.0.5
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1 LGG-04-1   6.576402               8
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1 LGG-04-1   6.693011               2
#                                     seurat_clusters harmony_clusters
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1              14               14
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1               0                0
#                                      copykat.pred
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1  not.defined
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1      diploid

obj_copykat <- RunUMAP(obj_copykat, reduction = "rpca", dims = 1:30)
# obj_copykat <- RunUMAP(obj_copykat, reduction = "pca", dims = 1:30)

p = DimPlot(obj_copykat, 
        reduction = "umap", 
        group.by = "copykat.pred", 
        label = TRUE, pt.size = 0.8) +
        ggtitle("copykat Predicted Aneuploidy Status")

ggsave(filename = file.path(output_dir, "copykat_UMAP_predicted_aneuploidy.png"),
       plot = p, 
       width = 7, 
       height = 6, 
       dpi = 300)

#####################################################################
#####################################################################

# str(copykat_result, max.level = 1)
# List of 3
# $ prediction :'data.frame':  499 obs. of  2 variables:
# $ CNAmat     :'data.frame':  12167 obs. of  472 variables:
# $ hclustering:List of 7
#  ..- attr(*, "class")= chr "hclust"

# str(copykat_result$prediction, max.level = 1)
#'data.frame': 499 obs. of  2 variables:
# $ cell.names  : chr  "LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1" "LGG-04-1_LGG-04-1_AAGTCGTCAGGCCCTA-1" "LGG-04-1_LGG-04-1_AATAGAGAGACCTCAT-1" "LGG-04-1_LGG-04-1_AATGACCTCAGTCTTT-1" ...
# $ copykat.pred: chr  "diploid" "aneuploid" "aneuploid" "diploid" ...
 
# str(copykat_result$CNAmat, max.level = 1)
#'data.frame': 12167 obs. of  472 variables:
# $ chrom                               : int  1 1 1 1 1 1 1 1 1 1 ...
# $ chrompos                            : num  1042457 1265484 1519859 1826619 2058465 ...
# $ abspos                              : num  1042457 1265484 1519859 1826619 2058465 ...

# str(copykat_result$hclustering, max.level = 1)
#ist of 7
# $ merge      : int [1:468, 1:2] -56 -302 -314 -244 -32 -7 -430 -251 -460 -9 ...
# $ height     : num [1:468] 3.4 3.54 3.64 3.68 3.71 ...
# $ order      : int [1:469] 16 10 90 1 207 18 454 43 123 115 ...
# $ labels     : chr [1:469] "LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1" "LGG.04.1_LGG.04.1_AAGTCGTCAGGCCCTA.1" "LGG.04.1_LGG.04.1_AATAGAGAGACCTCAT.1" "LGG.04.1_LGG.04.1_AATGACCTCAGTCTTT.1" ...
# $ method     : chr "ward.D"
# $ call       : language hclust(d = parallelDist::parDist(t(mat.adj), threads = n.cores, method = distance),      method = "ward.D")
# $ dist.method: chr "euclidean"
# - attr(*, "class")= chr "hclust"

#####################################################################
#####################################################################

# Check structure

dim(copykat_result$CNAmat)
head(rownames(copykat_result$CNAmat))   # gene names
head(colnames(copykat_result$CNAmat))   # cell barcodes
print(copykat_result$CNAmat[1:5, 1:5])

# CNA = copykat_result$CNAmat
# dim(CNA)
# to keep only the cells that are called "diploid" or "aneuploid" for subsequent analysis

# head(colnames(CNA), 6)
# [1] "chrom"                               
# [2] "chrompos"                            
# [3] "abspos"                              
# [4] "LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1"
# [5] "LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1"
# [6] "LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1"

#####################################################################
##################################################################### the output files are the following : 

# 📄 RPCA_copykat_prediction.txt
# What it is: The main prediction table.
# Content: Each row corresponds to a single cell and includes:
# Cell name
# copykat.pred: "diploid", "aneuploid", or "not.defined"
# CNV cluster assignment
# Possibly confidence scores
# You can load it with:
# pred_df <- read.table("RPCA_copykat_prediction.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 🧬 RPCA_copykat_CNA_results.txt
# What it is: Processed Copy Number Alteration (CNA) results, per cell and genomic region.
# Content: This is a matrix or table of CNV calls across the genome per cell.

# 🧬 RPCA_copykat_CNA_results.seg
# What it is: Segmented CNV results in SEG format, suitable for IGV or GISTIC2.
# Fields usually include: Chromosome, start, end, cell name, and segment mean log2 ratio.

# 🧬 RPCA_copykat_CNA_raw_results_gene_by_cell.txt
# What it is: The raw gene-by-cell CNV matrix, before smoothing.
# Content: A matrix with:
# Rows = genes
# Columns = cells
# Values = inferred log2 CNV

#####################################################################
#####################################################################

# The first 3 columns in the CNV matrix are the genomic coordinates. Rows are 220KB bins in genomic orders. 
# It also output CNVs indexed by gene names and other information such as chromosome names, start and end positions, 
# G staining band etc

#  chrom chrompos  abspos LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1
#1     1  1042457 1042457                           0.04784874
#2     1  1265484 1265484                           0.04784874
#3     1  1519859 1519859                           0.04784874
#4     1  1826619 1826619                           0.04784874
#5     1  2058465 2058465                           0.04784874
#6     1  2280372 2280372                           0.04784874
#  LGG.04.1_LGG.04.1_AAGTCGTCAGGCCCTA.1
#1                           -0.0824941
#2                           -0.0824941
#3                           -0.0824941
#4                           -0.0824941
#5                           -0.0824941
#6                           -0.0824941

# chrompos → for plotting within individual chromosomes.
# abspos → for genome-wide plots (e.g. heatmaps, Manhattan plots) where x-axis is continuous across chromosomes.

#########################################################
#########################################################

# 🔍 Visualization Effect:
# Range Description Bins
# -1 to -0.4  Strong deletions  50
# -0.4 to -0.2  Mild deletions  150
# -0.2 to 0.2 Near-diploid (normal variation region)  600
# 0.2 to 0.4  Mild amplifications 150
# 0.4 to 1  Strong amplifications 50

# 📌 Summary of Typical Thresholds:
# Platform  Gain Threshold  Loss Threshold
# SNP array / WGS log2 > 0.3–0.4  log2 < -0.3–0.4
# Absolute CN CN > 2  CN < 2
# Single-cell RNA log2 > 0.2–0.3  log2 < -0.2–0.3

#########################################################
#########################################################
#########################################################
######################################################### output files

# PREDICTIONS

# prediction_df <- read.table("RPCA_copykat_prediction.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# head(prediction_df)
#                            cell.names copykat.pred
# 1 LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1    aneuploid
# 2 LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1    aneuploid
# 3 LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1    aneuploid
# 4 LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1      diploid
# 5 LGG-04-1_LGG-04-1_AAGTCGTCAGGCCCTA-1      diploid
# 6 LGG-04-1_LGG-04-1_AATAGAGAGACCTCAT-1      diploid
 
# SEGMENTATION

# seg_df <- read.table("RPCA_copykat_CNA_results.seg", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# head(seg_df)
#                                    ID chrom loc.start  loc.end num.mark
# 1 LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1     1   1042457  3230151       10
# 2 LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1     1   3437000  3645603        2
# 3 LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1     1   4091348 10808784       31
# 4 LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1     1  11040607 17777104       23
# 5 LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1     1  17987772 23563511       26
# 6 LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1     1  24139449 27685740       16
#       seg.mean
# 1 -0.1154144144
# 2 -0.0688295141
# 3 -0.0002964707
# 4  0.0390021256
# 5  0.1289415969
# 6  0.1389694081

# RESULTS

# cna_results_df <- read.table("RPCA_copykat_CNA_results.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# head(cna_results_df)

# Show first 4 rows and 5 columns
# head(cna_results_df[, 1:4], 5)
#  chrom chrompos  abspos LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1
# 1     1  1042457 1042457                           -0.1154144
# 2     1  1265484 1265484                           -0.1154144
# 3     1  1519859 1519859                           -0.1154144
# 4     1  1826619 1826619                           -0.1154144
#  LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1 LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1
# 1                          -0.05418733                           0.02951918
# 2                          -0.05418733                           0.02951918
# 3                          -0.05418733                           0.02951918
# 4                          -0.05418733                           0.02951918

# Column Name Description
# chrom Chromosome number where the data point (usually a gene or genomic bin) is located.
# chrompos  Chromosome-relative genomic coordinate (start or mid-point of gene/bin).
# abspos  Absolute genomic position across the whole genome (used for continuous plotting). May be the same as chrompos unless chromosome offsets were added.
# Columns 4+  Each of the remaining columns corresponds to one single cell, with a value representing its log2 copy number ratio at that genomic position.

# 🧠 Interpretation of Cell Values:
# Each number is a log2 ratio of inferred copy number:

# Values ≈ 0 → likely diploid (normal)
# Values > 0.2 → likely gain
# Values < -0.2 → likely loss

# The files copykat_result$CNAmat and cna_results_df are identical.

# print(copykat_result$CNAmat[1:5, 1:5])
#  chrom chrompos  abspos LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1
# 1     1  1042457 1042457                           -0.1154144
# 2     1  1265484 1265484                           -0.1154144
# 3     1  1519859 1519859                           -0.1154144
# 4     1  1826619 1826619                           -0.1154144
# 5     1  2058465 2058465                           -0.1154144
#  LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1
# 1                          -0.05418733
# 2                          -0.05418733
# 3                          -0.05418733
# 4                          -0.05418733
# 5                          -0.05418733

# head(cna_results_df[1:5, 1:5])
#  chrom chrompos  abspos LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1
# 1     1  1042457 1042457                           -0.1154144
# 2     1  1265484 1265484                           -0.1154144
# 3     1  1519859 1519859                           -0.1154144
# 4     1  1826619 1826619                           -0.1154144
# 5     1  2058465 2058465                           -0.1154144
#  LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1
# 1                          -0.05418733
# 2                          -0.05418733
# 3                          -0.05418733
# 4                          -0.05418733
# 5                          -0.05418733

# Read the Gene file (cr, start, end)

# cna_raw_df <- read.table("RPCA_copykat_CNA_raw_results_gene_by_cell.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# head(colnames(cna_raw_df), 10)
# [1] "abspos"                              
# [2] "chromosome_name"                     
# [3] "start_position"                      
# [4] "end_position"                        
# [5] "ensembl_gene_id"                     
# [6] "hgnc_symbol"                         
# [7] "band"                                
# [8] "LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1"

# display only first 10 cells for readability
#         chromosome_name start_position end_position ensembl_gene_id hgnc_symbol

# head(cna_raw_df[, 1:10]) 
#   abspos chromosome_name start_position end_position ensembl_gene_id
# 1  951756               1         944204       959309 ENSG00000188976
# 2  999567               1         998962      1000172 ENSG00000188290
# 3 1007839               1        1001138      1014541 ENSG00000187608
# 4 1224469               1        1216908      1232031 ENSG00000078808
# 5 1233653               1        1232265      1235041 ENSG00000176022
# 6 1263897               1        1253909      1273885 ENSG00000160087
#  hgnc_symbol   band LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1
# 1       NOC2L p36.33                           -0.1315771
# 2        HES4 p36.33                           -0.1315771
# 3       ISG15 p36.33                           -0.1315771
# 4        SDF4 p36.33                           -0.1315771
# 5     B3GALT6 p36.33                           -0.1315771
# 6      UBE2J2 p36.33                           -0.1315771
#  LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1 LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1
# 1                          -0.07579418                          0.002886359
# 2                          -0.07579418                          0.002886359
# 3                          -0.07579418                          0.002886359
# 4                          -0.07579418                          0.002886359
# 5                          -0.07579418                          0.002886359
# 6                          -0.07579418                          0.002886359

# 📦 How CopyKAT generates the .seg file:
# 🔧 Internally:
# Gene expression matrix is smoothed over genomic bins or windows (e.g., win.size = 25).
# CopyKAT performs sliding window averaging across genes ordered by genomic position.
# It infers relative log2 copy number values per window per cell.
# For each cell, CopyKAT segments contiguous genomic regions with similar CNV signal into:
# Chromosome
# Start and end position
# Segment mean (log2 CNV estimate)
# The output is formatted in SEG format, which is widely used for genome visualization tools (like IGV or GISTIC).

# 📄 Format of the .seg file
# A typical .seg file has these columns:

# Column Name Description
# ID  Cell barcode
# chrom Chromosome (e.g., 1, 2, X, Y)
# loc.start Start position of the segment (genomic coordinate)
# loc.end End position of the segment
# num.mark  Number of genes/markers in the segment
# seg.mean  Mean log2 CNV ratio for the segment

# to keep only the cells that are assigned to diploid or aneuploid phenotype
# also copykat_result$CNAmat has column names with dots (.) instead of hyphens (-), 
# and also contains extra columns like "chrom", "chrompos", and "abspos" that are not cells.

# head(obj_copykat@meta.data, 4)
#                                     orig.ident nCount_RNA nFeature_RNA
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1   LGG-04-1        517          265
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1   LGG-04-1       7999         2722
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1   LGG-04-1       9911         3005
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1   LGG-04-1      16196         3761
#                                       sample percent.mt RNA_snn_res.0.5
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1 LGG-04-1   6.576402               8
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1 LGG-04-1   8.738592               0
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1 LGG-04-1   8.858844               0
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1 LGG-04-1   6.693011               2
#                                     seurat_clusters harmony_clusters
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1              14               14
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1               0                0
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1               0                0
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1               0                0
#                                     copykat.pred
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1  not.defined
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1    aneuploid
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1    aneuploid
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1    aneuploid

# head(colnames(copykat_result$CNAmat), 6)
#[1] "chrom"                               
#[2] "chrompos"                            
#[3] "abspos"                              
#[4] "LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1"
#[5] "LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1"
#[6] "LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1"

##################################################################### !
##################################################################### !
# copyKat introduces . in the name of the cells (to correct it ...)
##################################################################### !
##################################################################### !

# head(colnames(copykat_result$prediction), 6)
# [1] "cell.names"   "copykat.pred"
# head((copykat_result$prediction), 6)
#                                                               cell.names
# LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1 LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1
# LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1 LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1
# LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1 LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1
# LGG.04.1_LGG.04.1_AAGTACCTCTCAACCC.1 LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1
# LGG.04.1_LGG.04.1_AAGTCGTCAGGCCCTA.1 LGG-04-1_LGG-04-1_AAGTCGTCAGGCCCTA-1
# LGG.04.1_LGG.04.1_AATAGAGAGACCTCAT.1 LGG-04-1_LGG-04-1_AATAGAGAGACCTCAT-1
#                                     copykat.pred
# LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1    aneuploid
# LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1    aneuploid
# LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1    aneuploid
# LGG.04.1_LGG.04.1_AAGTACCTCTCAACCC.1      diploid
# LGG.04.1_LGG.04.1_AAGTCGTCAGGCCCTA.1      diploid
# LGG.04.1_LGG.04.1_AATAGAGAGACCTCAT.1      diploid

# head(copykat_pred_filtered, 3)
#                                                              cell.names
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1 LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1 LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1 LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1
#                                     copykat.pred
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1    aneuploid
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1    aneuploid
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1    aneuploid

# dim(copykat_pred)
# [1] 1000    2
# dim(copykat_pred_filtered)
# [1] 947   2

# head(colnames(copykat_result$CNAmat))
# [1] "chrom"                               
# [2] "chrompos"                            
# [3] "abspos"                              
# [4] "LGG.04.1_LGG.04.1_AACGTCACAACGGCCT.1"
# [5] "LGG.04.1_LGG.04.1_AAGCATCTCCTCGATC.1"
# [6] "LGG.04.1_LGG.04.1_AAGGAATCACGCCACA.1"
# head(rownames(copykat_pred_filtered))
# [1] "LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1"
# [2] "LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1"
# [3] "LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1"
# [4] "LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1"
# [5] "LGG-04-1_LGG-04-1_AAGTCGTCAGGCCCTA-1"
# [6] "LGG-04-1_LGG-04-1_AATAGAGAGACCTCAT-1"

# copykat_pred <- copykat_result$prediction
# copykat_pred_filtered <- subset(copykat_pred, copykat.pred %in% c("aneuploid", "diploid"))
# dim(copykat_pred_filtered)

# copykat_pred_filtered$cell.names.dot <- gsub("-", ".", copykat_pred_filtered$cell.names)
# rownames(copykat_pred_filtered) <- copykat_pred_filtered$cell.names.dot
# rownames(copykat_pred_filtered) <- copykat_pred_filtered$cell.names

# CNAmatrix_filtered <- copykat_result$CNAmat[, c("chrom", "chrompos", "abspos", rownames(copykat_pred_filtered))]

# Step 3: Check dimensions
# dim(CNAmatrix_filtered)
# [1] 12167   950
# dim(copykat_result$CNAmat)
# [1] 12167   950

#####################################################################
#####################################################################
#########################################################
######################################################### for simplicity : the variable is CNA
#####################################################################
#####################################################################

# CNA = CNAmatrix_filtered

# Count cells in diploid vs aneuploid groups
# aneuploid_count <- sum(copykat_result$prediction$copykat.pred == "aneuploid")
# diploid_count   <- sum(copykat_result$prediction$copykat.pred== "diploid")
# total_cells     <- length(copykat_result$prediction$copykat.pred)

# logR_mat <- CNA[, 4:ncol(CNA)]
# summary(logR_mat)

# Print the summary on one file
# logR_all <- as.vector(as.matrix(logR_mat))
# logR_all <- logR_all[!is.na(logR_all)]

# Define output file path
# summary_file <- file.path(output_dir, "copykat_logR_global_summary.txt")

# Write everything to file
# sink(summary_file)
# cat("Summary of all log2R values across all cells:\n\n")
# print(summary(logR_all))

# cat("\n\nCell classification counts:\n")
# cat("  Aneuploid cells: ", aneuploid_count, "\n")
# cat("  Diploid cells:   ", diploid_count, "\n")
# cat("  Total cells:     ", total_cells, "\n")
# sink()

# # Print summary
# summary(logR_all)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.361213 -0.037438 -0.000065  0.000000  0.034820  0.419477 

##################################################################### COMPLEXHEATMAP 1
##################################################################### 

----------------------------
# Step 1: Filter predictions to diploid and aneuploid
----------------------------
copykat_pred <- copykat_result$prediction
copykat_pred_filtered <- subset(copykat_pred, copykat.pred %in% c("aneuploid", "diploid"))
copykat_pred_filtered$cell.names.dot <- gsub("-", ".", copykat_pred_filtered$cell.names)
rownames(copykat_pred_filtered) <- copykat_pred_filtered$cell.names.dot

----------------------------
# Step 2: Subset CNA matrix
----------------------------
CNAmatrix_filtered <- copykat_result$CNAmat[, c("chrom", "chrompos", "abspos", rownames(copykat_pred_filtered))]
CNA <- CNAmatrix_filtered  # working variable

----------------------------
# Step 3: Count classification stats
----------------------------
aneuploid_count <- sum(copykat_pred_filtered$copykat.pred == "aneuploid")
diploid_count   <- sum(copykat_pred_filtered$copykat.pred == "diploid")
total_cells     <- ncol(CNA) - 3  # exclude metadata cols

----------------------------
# Step 4: LogR summary
----------------------------
logR_mat <- CNA[, 4:ncol(CNA)]
logR_all <- as.vector(as.matrix(logR_mat))
logR_all <- logR_all[!is.na(logR_all)]

# Save summary stats to file
summary_file <- file.path(output_dir, "copykat_logR_global_summary.txt")
sink(summary_file)
cat("Summary of all log2R values across all cells:\n\n")
print(summary(logR_all))
cat("\n\nCell classification counts:\n")
cat("  Aneuploid cells: ", aneuploid_count, "\n")
cat("  Diploid cells:   ", diploid_count, "\n")
cat("  Total cells:     ", total_cells, "\n")
sink()

----------------------------
# Step 5: Heatmap preparation
----------------------------
cna_matrix <- CNA[, 4:ncol(CNA)]

# Align prediction vector
pred_vector <- copykat_pred_filtered$copykat.pred
names(pred_vector) <- rownames(copykat_pred_filtered)
pred_vector <- pred_vector[colnames(cna_matrix)]
stopifnot(length(pred_vector) == ncol(cna_matrix))

# RowSideColors: Prediction
rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
pred_colors <- rbPal5(2)[as.numeric(factor(pred_vector))]
cells <- matrix(pred_colors, nrow = 1)

# ColSideColors: Chromosome annotation
chromosomes <- as.character(CNA$chrom)
chr_levels <- c(as.character(1:22), "X", "Y")
chromosomes <- factor(chromosomes, levels = chr_levels)
chr <- as.numeric(chromosomes) %% 2 + 1
rbPal1 <- colorRampPalette(c('black', 'grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR, CHR)

----------------------------
# Step 6: Plot and save heatmap
----------------------------

# Heatmap color breaks
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
col_breaks <- c(
  seq(-1, -0.4, length = 50),
  seq(-0.4, -0.2, length = 150),
  seq(-0.2, 0.2, length = 600),
  seq(0.2, 0.4, length = 150),
  seq(0.4, 1, length = 50)
)

heatmap_file <- file.path(output_dir, "copykat_genome_wide_complexheatmap1.png")
png(heatmap_file, width = 2400, height = 1800, res = 300)

heatmap.3(
  t(cna_matrix),
  dendrogram = "r",
  distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"),
  hclustfun = function(x) hclust(x, method = "ward.D2"),
  ColSideColors = chr1,
  RowSideColors = cells,
  Colv = NA,
  Rowv = TRUE,
  labRow = NA,
  labCol = NA,
  notecol = "black",
  col = my_palette,
  breaks = col_breaks,
  key = TRUE,
  keysize = 1,
  density.info = "none",
  trace = "none",
  cexRow = 0.1,
  cexCol = 0.1,
  cex.main = 1,
  cex.lab = 0.1,
  symm = FALSE,
  symkey = FALSE,
  symbreaks = TRUE,
  cex = 1,
  margins = c(10, 10)
)

# Add legend for predictions
legend(
  "topright",
  legend = paste("pred.", names(table(pred_vector)), sep = ""),
  pch = 15,
  col = unique(pred_colors),
  cex = 0.6,
  bty = "n"
)

dev.off()

#####################################################################
#####################################################################
#####################################################################
##################################################################### COMPLEXHEATMAP 2

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(parallelDist)
library(grid)

# Create output directory if needed
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

--
# Filter predictions and subset CNA
--
copykat_pred <- copykat_result$prediction
copykat_pred_filtered <- subset(copykat_pred, copykat.pred %in% c("aneuploid", "diploid"))
copykat_pred_filtered$cell.names.dot <- gsub("-", ".", copykat_pred_filtered$cell.names)
rownames(copykat_pred_filtered) <- copykat_pred_filtered$cell.names.dot

# Subset CNA matrix
CNAmatrix_filtered <- copykat_result$CNAmat[, c("chrom", "chrompos", "abspos", rownames(copykat_pred_filtered))]
CNA <- CNAmatrix_filtered
cna_matrix <- CNA[, 4:ncol(CNA)]

# Align prediction vector
pred_vector <- copykat_pred_filtered$copykat.pred
names(pred_vector) <- rownames(copykat_pred_filtered)
pred_vector <- pred_vector[colnames(cna_matrix)]
stopifnot(length(pred_vector) == ncol(cna_matrix))

--
# Custom sharp color palette
--
sharp_colors <- c(
  "#0000FF", "#6699FF", "#CCDDFF", "#FFFFFF", "#FFDDCC", "#FF9966", "#FF0000"
)
col_fun <- colorRamp2(
  breaks = c(-1, -0.4, -0.2, 0, 0.2, 0.4, 1),
  colors = sharp_colors
)

--
# Chromosome metadata
--
chromosomes <- as.character(CNA$chrom)
chr_levels <- c(as.character(1:22), "X", "Y")
chromosomes <- factor(chromosomes, levels = chr_levels)
chr_levels <- chr_levels[chr_levels %in% unique(chromosomes)]

--
# Prediction annotation bar
--
pred_color_map <- brewer.pal(n = 8, name = "Dark2")[2:1]
names(pred_color_map) <- unique(pred_vector)

anno <- rowAnnotation(
  Prediction = pred_vector,
  col = list(Prediction = pred_color_map),
  show_annotation_name = TRUE,
  annotation_name_rot = 0,
  show_legend = FALSE
)

--
# Row clustering
--
row_dend <- hclust(parallelDist::parDist(t(cna_matrix), threads = 4, method = "euclidean"),
                   method = "ward.D2")

--
# Generate heatmaps per chromosome
--
ht_list <- list()
for (chr in chr_levels) {
  chr_indices <- which(chromosomes == chr)
  if (length(chr_indices) > 0) {
    chr_data <- t(cna_matrix)[, chr_indices, drop = FALSE]
    ht_list[[chr]] <- Heatmap(
      chr_data,
      name = paste0("chr", chr),
      col = col_fun,
      cluster_rows = row_dend,
      cluster_columns = FALSE,
      column_title = paste0("chr", chr),
      column_title_rot = 90,
      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
      show_row_names = FALSE,
      show_column_names = FALSE,
      show_heatmap_legend = FALSE
    )
  }
}

# Combine all heatmaps
ht_combined <- Reduce(`+`, ht_list) + anno

--
# Define legends
--
lgd_list <- list(
  Legend(
    title = "Copy Number",
    at = c(-1, -0.4, -0.2, 0, 0.2, 0.4, 1),
    labels = c("Loss", "Mod. Loss", "Slight Loss", "Neutral", "Slight Gain", "Mod. Gain", "Gain"),
    col_fun = col_fun
  ),
  Legend(
    title = "Prediction",
    labels = names(pred_color_map),
    legend_gp = gpar(fill = pred_color_map)
  )
)

--
# Save the figure
--
heatmap_file2 <- file.path(output_dir, "copykat_genome_wide_complexheatmap2.png")
png(heatmap_file2, width = 3600, height = 1800, res = 300)

draw(ht_combined, 
     padding = unit(c(15, 2, 2, 2), "mm"), 
     heatmap_legend_list = lgd_list, 
     heatmap_legend_side = "right")

dev.off()

#####################################################################
#####################################################################
#####################################################################
##################################################################### COMPLEXHEATMAPS : multiple methods
######################################################### logR : -1 to +1
######################################################### separation for moderate vs. deep CNV
#########################################################
######################################################### MULTIPLE clustering methods METHODS

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(parallelDist)
library(proxy)
library(grid)

# Ensure output_dir exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

---
# Filter predictions and subset CNA matrix
---
copykat_pred <- copykat_result$prediction
copykat_pred_filtered <- subset(copykat_pred, copykat.pred %in% c("aneuploid", "diploid"))
copykat_pred_filtered$cell.names.dot <- gsub("-", ".", copykat_pred_filtered$cell.names)
rownames(copykat_pred_filtered) <- copykat_pred_filtered$cell.names.dot

# Subset CNA matrix to filtered cells
CNAmatrix_filtered <- copykat_result$CNAmat[, c("chrom", "chrompos", "abspos", rownames(copykat_pred_filtered))]
CNA <- CNAmatrix_filtered
cna_matrix <- CNA[, 4:ncol(CNA)]

# Align prediction vector
pred_vector <- copykat_pred_filtered$copykat.pred
names(pred_vector) <- rownames(copykat_pred_filtered)
pred_vector <- pred_vector[colnames(cna_matrix)]
stopifnot(length(pred_vector) == ncol(cna_matrix))

---
# Define color palette and metadata
---
sharp_colors <- c("#0000FF", "#6699FF", "#CCDDFF", "#FFFFFF", "#FFDDCC", "#FF9966", "#FF0000")
col_fun <- colorRamp2(
  breaks = c(-1, -0.4, -0.2, 0, 0.2, 0.4, 1),
  colors = sharp_colors
)

# Chromosome information
chromosomes <- as.character(CNA$chrom)
chr_levels <- c(as.character(1:22), "X", "Y")
chr_levels <- chr_levels[chr_levels %in% unique(chromosomes)]
chr_factor <- factor(chromosomes, levels = chr_levels)

# Prediction color map
pred_colors <- brewer.pal(n = 8, name = "Dark2")[2:1]
names(pred_colors) <- unique(pred_vector)

# Transpose matrix: rows = cells
data_matrix <- t(cna_matrix)

---
# Iterate over distance methods
---
methods <- c("euclidean", "manhattan", "canberra", "minkowski", "cosine")

for (method in methods) {
  cat("Processing:", method, "\n")
  
  # Compute row distances (i.e., for cells)
  if (method == "cosine") {
    d <- proxy::dist(data_matrix, method = "cosine")
  } else {
    d <- parallelDist::parDist(data_matrix, method = method, threads = 4)
  }
  
  row_dend <- hclust(d, method = "ward.D2")
  
  # Build chromosome-wise heatmaps
  ht_list <- list()
  for (chr in chr_levels) {
    chr_indices <- which(chromosomes == chr)
    if (length(chr_indices) > 0) {
      chr_data <- data_matrix[, chr_indices, drop = FALSE]
      ht_list[[chr]] <- Heatmap(
        chr_data,
        name = paste0("chr", chr),
        col = col_fun,
        cluster_rows = row_dend,
        cluster_columns = FALSE,
        column_title = paste0("chr", chr),
        column_title_rot = 90,
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE
      )
    }
  }

  # Row annotation for predictions
  anno <- rowAnnotation(
    Prediction = pred_vector,
    col = list(Prediction = pred_colors),
    show_annotation_name = TRUE,
    annotation_name_rot = 0,
    show_legend = FALSE
  )

  # Combine heatmaps
  ht_combined <- Reduce(`+`, ht_list) + anno

  # Legend list
  lgd_list <- list(
    Legend(
      title = "Copy Number", 
      at = c(-1, -0.4, -0.2, 0, 0.2, 0.4, 1),
      labels = c("Loss (-1)", "Mod. Loss (-0.4)", "Slight Loss (-0.2)",
                 "Neutral (0)", "Slight Gain (0.2)", "Mod. Gain (0.4)", "Gain (1)"),
      col_fun = col_fun
    ),
    Legend(
      title = "Prediction",
      labels = names(pred_colors),
      legend_gp = gpar(fill = pred_colors)
    )
  )

  # Save output
  out_file <- file.path(output_dir, paste0("copykat_genome_wide_complexheatmap2_", method, ".png"))
  png(out_file, width = 2600, height = 1800, res = 300)

  draw(ht_combined, padding = unit(c(15, 2, 2, 2), "mm"), 
       heatmap_legend_list = lgd_list,
       heatmap_legend_side = "right")

  dev.off()
}

#####################################################################
#####################################################################
######################################################## logR : -0.5 to 0.5
######################################################## multiple methods
#####################################################################
#####################################################################

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(parallelDist)
library(proxy)
library(RColorBrewer)
library(grid)

------
# Step 1: Filter predictions & subset CNA
------
copykat_pred <- copykat_result$prediction
copykat_pred_filtered <- subset(copykat_pred, copykat.pred %in% c("aneuploid", "diploid"))
copykat_pred_filtered$cell.names.dot <- gsub("-", ".", copykat_pred_filtered$cell.names)
rownames(copykat_pred_filtered) <- copykat_pred_filtered$cell.names.dot

# Subset CNA matrix
CNAmatrix_filtered <- copykat_result$CNAmat[, c("chrom", "chrompos", "abspos", rownames(copykat_pred_filtered))]
CNA <- CNAmatrix_filtered
cna_matrix <- CNA[, 4:ncol(CNA)]

# Align prediction vector
pred_vector <- copykat_pred_filtered$copykat.pred
names(pred_vector) <- rownames(copykat_pred_filtered)
pred_vector <- pred_vector[colnames(cna_matrix)]
stopifnot(length(pred_vector) == ncol(cna_matrix))

------
# Step 2: Define color scales and metadata
------
col_fun <- colorRamp2(
  breaks = c(-0.5, -0.25, 0, 0.2, 0.5),
  colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#FDAE6B", "#E6550D")
)

# Chromosome info
chromosomes <- as.character(CNA$chrom)
chr_levels <- c(as.character(1:22), "X", "Y")
chr_levels <- chr_levels[chr_levels %in% unique(chromosomes)]

# Prediction annotation colors
pred_colors <- brewer.pal(n = 8, name = "Dark2")[2:1]
names(pred_colors) <- unique(pred_vector)

# Transpose CNA matrix: rows = cells
data_matrix <- t(cna_matrix)

------
# Step 3: Iterate over clustering methods
------
methods <- c("euclidean", "manhattan", "canberra", "minkowski", "cosine")

for (method in methods) {
  cat("Processing:", method, "\n")
  
  # Compute distance matrix for cells
  if (method == "cosine") {
    d <- proxy::dist(data_matrix, method = "cosine")
  } else {
    d <- parallelDist::parDist(data_matrix, method = method, threads = 4)
  }
  
  row_dend <- hclust(d, method = "ward.D2")
  
  # Create heatmap per chromosome
  ht_list <- list()
  for (chr in chr_levels) {
    chr_indices <- which(chromosomes == chr)
    if (length(chr_indices) > 0) {
      chr_data <- data_matrix[, chr_indices, drop = FALSE]
      
      ht_list[[chr]] <- Heatmap(
        chr_data,
        name = paste0("chr", chr),
        col = col_fun,
        cluster_rows = row_dend,
        cluster_columns = FALSE,
        column_title = paste0("chr", chr),
        column_title_rot = 90,
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE
      )
    }
  }
  
  # Row annotation (predicted CopyKAT status)
  anno <- rowAnnotation(
    Prediction = pred_vector,
    col = list(Prediction = pred_colors),
    show_annotation_name = TRUE,
    annotation_name_rot = 0,
    show_legend = FALSE
  )
  
  # Combine all chromosome heatmaps
  ht_combined <- Reduce(`+`, ht_list) + anno
  
  # Define legends
  lgd_list <- list(
    Legend(
      title = "Copy Number", 
      at = c(-0.5, -0.25, 0, 0.2, 0.5),
      labels = c("Loss (-0.5)", "Slight Loss (-0.25)", "Neutral (0)", 
                 "Slight Gain (0.2)", "Gain (0.5)"),
      col_fun = col_fun
    ),
    Legend(
      title = "Prediction", 
      labels = names(pred_colors),
      legend_gp = gpar(fill = pred_colors)
    )
  )
  
  # Save heatmap
  heatmap_file <- file.path(output_dir, paste0("copykat_genome_wide_complexheatmap3_", method, ".png"))
  png(heatmap_file, width = 2600, height = 1800, res = 300)
  
  draw(
    ht_combined,
    padding = unit(c(15, 2, 2, 2), "mm"),
    heatmap_legend_list = lgd_list,
    heatmap_legend_side = "right"
  )
  
  dev.off()
}

#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
# Chromosome-wise Heatmaps with Updated Palette

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(parallelDist)
library(proxy)
library(RColorBrewer)
library(grid)


------
# Step 1: Filter predictions & subset CNA
------
copykat_pred <- copykat_result$prediction
copykat_pred_filtered <- subset(copykat_pred, copykat.pred %in% c("aneuploid", "diploid"))
copykat_pred_filtered$cell.names.dot <- gsub("-", ".", copykat_pred_filtered$cell.names)
rownames(copykat_pred_filtered) <- copykat_pred_filtered$cell.names.dot

# Subset CNA matrix
CNAmatrix_filtered <- copykat_result$CNAmat[, c("chrom", "chrompos", "abspos", rownames(copykat_pred_filtered))]
CNA <- CNAmatrix_filtered
cna_matrix <- CNA[, 4:ncol(CNA)]

# Align prediction vector
pred_vector <- copykat_pred_filtered$copykat.pred
names(pred_vector) <- rownames(copykat_pred_filtered)
pred_vector <- pred_vector[colnames(cna_matrix)]
stopifnot(length(pred_vector) == ncol(cna_matrix))

------
# Step 2: Define color scales and metadata
------
col_fun <- colorRamp2(
  breaks = c(-0.5, -0.25, 0, 0.2, 0.5),
  colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#FDAE6B", "#E6550D")
)

# Chromosome info
chromosomes <- as.character(CNA$chrom)
chr_levels <- c(as.character(1:22), "X", "Y")
chr_levels <- chr_levels[chr_levels %in% unique(chromosomes)]

# Prediction annotation colors
pred_colors <- brewer.pal(n = 8, name = "Dark2")[2:1]
names(pred_colors) <- unique(pred_vector)

# Transpose CNA matrix: rows = cells
data_matrix <- t(cna_matrix)

------
# Step 3: Iterate over clustering methods
------
methods <- c("euclidean", "manhattan", "canberra", "minkowski", "cosine")

for (method in methods) {
  cat("Processing:", method, "\n")
  
  # Compute distance matrix for cells
  if (method == "cosine") {
    d <- proxy::dist(data_matrix, method = "cosine")
  } else {
    d <- parallelDist::parDist(data_matrix, method = method, threads = 4)
  }
  
  row_dend <- hclust(d, method = "ward.D2")
  
  # Create heatmap per chromosome
  ht_list <- list()
  for (chr in chr_levels) {
    chr_indices <- which(chromosomes == chr)
    if (length(chr_indices) > 0) {
      chr_data <- data_matrix[, chr_indices, drop = FALSE]
      
      ht_list[[chr]] <- Heatmap(
        chr_data,
        name = paste0("chr", chr),
        col = col_fun,
        cluster_rows = row_dend,
        cluster_columns = FALSE,
        column_title = paste0("chr", chr),
        column_title_rot = 90,
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE
      )
    }
  }
  
  # Row annotation (predicted CopyKAT status)
  anno <- rowAnnotation(
    Prediction = pred_vector,
    col = list(Prediction = pred_colors),
    show_annotation_name = TRUE,
    annotation_name_rot = 0,
    show_legend = FALSE
  )
  
  # Combine all chromosome heatmaps
  ht_combined <- Reduce(`+`, ht_list) + anno
  
  # Define legends
  lgd_list <- list(
    Legend(
      title = "Copy Number", 
      at = c(-0.5, -0.25, 0, 0.2, 0.5),
      labels = c("Loss (-0.5)", "Slight Loss (-0.25)", "Neutral (0)", 
                 "Slight Gain (0.2)", "Gain (0.5)"),
      col_fun = col_fun
    ),
    Legend(
      title = "Prediction", 
      labels = names(pred_colors),
      legend_gp = gpar(fill = pred_colors)
    )
  )
  
  # Save heatmap
  heatmap_file <- file.path(output_dir, paste0("copykat_genome_wide_complexheatmap4_", method, ".png"))
  png(heatmap_file, width = 2600, height = 1800, res = 300)
  
  draw(
    ht_combined,
    padding = unit(c(15, 2, 2, 2), "mm"),
    heatmap_legend_list = lgd_list,
    heatmap_legend_side = "right"
  )
  
  dev.off()
}

#####################################################################
#####################################################################
########################################################
######################################################## seg files 
########################################################
######################################################## CNV frequency
#####################################################################
#####################################################################


library(svpluscnv)
# library(Cairo)  

# Read the CopyKAT segmentation file
# seg <- read.table(paste0(illness,"_copykat_CNA_results.seg", sep=""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set expected column names
# colnames(seg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
# colnames(seg)
# [1] "ID"        "chrom"     "loc.start" "loc.end"   "num.mark"  "seg.mean"

# each cell is a sample :
# length(unique(seg$ID))
# [1] 947

prefix <- sub("_.*", "", illness)
print(prefix)  # Output: "LGG"

seg_file <- file.path(output_dir, paste0(prefix, "_copykat_CNA_results.seg"))

# Check if the file exists
if (file.exists(seg_file)) {
  message("The file exists: ", seg_file)
} else {
  message("The file does not exist: ", seg_file)
}

seg <- read.table(seg_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cnv = validate.cnv(seg)

# An object of class svcnvio from svpluscnv storing cnv data from 947 samples

cnv_freq <- cnv.freq(cnv, fc.pct = 0.2, ploidy = FALSE)

# str(cnv_freq)
# Formal class 'cnvfreq' [package "svpluscnv"] with 3 slots
#  ..@ freqsum:Classes ‘data.table’ and 'data.frame':  2936 obs. of  5 variables:
#  .. ..$ chr       : chr [1:2936] "chr1" "chr1" "chr1" "chr1" ...
#  .. ..$ start     : num [1:2936] 1042457 2042457 3042457 4042457 5042457 ...
#  .. ..$ end       : num [1:2936] 2042457 3042457 4042457 5042457 6042457 ...
#  .. ..$ freq.gains: num [1:2936] 0.00106 0.00106 0.02429 0 0 ...
#  .. ..$ freq.loss : num [1:2936] 0 0 0.0253 0 0 ...
#  .. ..- attr(*, ".internal.selfref")=<externalptr> 
#  ..@ bin.mat: num [1:2936, 1:947] 9.64e-05 9.64e-05 9.64e-05 9.64e-05 9.64e-05 ...
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:2936] "chr1_1042457_2042457" "chr1_2042457_3042457" "chr1_3042457_4042457" "chr1_4042457_5042457" ...
#  .. .. ..$ : chr [1:947] "LGG.03_LGG.03_AAAGAACGTAGCGTTT.1" "LGG.03_LGG.03_AAAGGATAGCTAGTTC.1" "LGG.03_LGG.03_AAAGGGCCAGCTTTCC.1" "LGG.03_LGG.03_AAAGTGATCAGCAGAG.1" ...
#  ..@ param  :List of 7

# cnv_freq@freqsum
#         chr     start       end  freq.gains  freq.loss
#      <char>     <num>     <num>       <num>      <num>
#   1:   chr1   1042457   2042457 0.001055966 0.00000000
#   2:   chr1   2042457   3042457 0.001055966 0.00000000
#   3:   chr1   3042457   4042457 0.024287223 0.02534319
#   4:   chr1   4042457   5042457 0.000000000 0.00000000
#   5:   chr1   5042457   6042457 0.000000000 0.00000000
#  ---                                                  
# 2932:  chr23 151517443 152517443 0.000000000 0.00000000
# 2933:  chr23 152517443 153517443 0.000000000 0.00000000
# 2934:  chr23 153517443 154517443 0.000000000 0.00000000
# 2935:  chr23 154517443 155517443 0.000000000 0.00000000
# 2936:  chr23 155517443 156517443 0.000000000 0.00000000

# cnv.freq(
#       cnv,
#       fc.pct = 0.2,
#       genome.v = "hg19",
#       ploidy = FALSE,
#       g.bin = 1,
#       sampleids = NULL,
#       cex.axis = 1,
#       cex.lab = 1,
#       label.line = -1.2,
#       plot = TRUE,
#       summary = TRUE,
#       verbose = TRUE

# Define the output file path
png_file <- file.path(output_dir, "copykat_genome_wide_cnv_frequency_cnvplussv.png")
# Save plot with transparency support to the specified folder
png(png_file, width = 1400, height = 600, res = 150)

# Run frequency calculation with hg38 genome and plotting
cnv.freq(
  cnv,                # validated CNV object
  fc.pct = 0.2,       # fold-change threshold for calling gain/loss
  genome.v = "hg38",  # ✅ specify hg38 instead of default hg19
  ploidy = FALSE,     # use absolute seg.mean thresholds
  g.bin = 1,          # genome binning resolution (Mb)
  sampleids = NULL,   # all samples
  cex.axis = 0.8,     # increase axis text size
  cex.lab = 0.8,      # increase label size
  label.line = -1.5,  # controls y-axis title offset
  plot = TRUE,        # enable default plotting
  summary = TRUE,
  verbose = TRUE
)

dev.off()

#####################################################################
#####################################################################
####################################### customized visualization
#####################################################################
#####################################################################

library(ggplot2)
library(dplyr)

# Extract and rename columns from freqsum slot
cnv_df <- cnv_freq@freqsum %>%
  dplyr::rename(
    chrom = chr,
    loc.start = start,
    loc.end = end,
    gain.freq = freq.gains,
    loss.freq = freq.loss
  )

# Compute midpoint in Mb
cnv_df$midpoint <- (cnv_df$loc.start + cnv_df$loc.end) / 2 / 1e6

# Clean and factor chromosome labels
cnv_df$chrom <- gsub("chr", "", cnv_df$chrom)
cnv_df$chrom <- factor(cnv_df$chrom, levels = as.character(1:23))  # Include chrX=23 if needed

# Improved plot
p = ggplot(cnv_df, aes(x = midpoint)) +
  geom_line(aes(y = gain.freq, color = "Gain"), size = 1.2) +
  geom_line(aes(y = -loss.freq, color = "Loss"), size = 1.2) +
  facet_wrap(~chrom, scales = "free_x", nrow = 4) +
  scale_color_manual(values = c("Gain" = "#E41A1C", "Loss" = "#377EB8")) +  # red for gain, blue for loss
  labs(
    title = "Copy Number Alteration Frequency (CopyKAT)",
    x = "Genomic Position (Mb)",
    y = "Gain (+) / Loss (-) Frequency",
    color = "CNA Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

# Save the ggplot to PDF in the specified folder
ggsave(
  filename = file.path(output_dir, "copykat_genome_wide_cnv_freq_plot_custom.pdf"),
  plot = p,
  width = 24,
  height = 12,
  bg = "transparent"
)

#####################################################################
#####################################################################

# Run genome alteration percentage
pct_change <- pct.genome.changed(cnv, fc.pct = 0.2)

# Save to text file (tab-delimited)
write.table(
  data.frame(pct_change),
  file = file.path(output_dir, "copykat_pct_genome_altered.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

#####################################################################
##################################################################### VIOLIN PLOT
######################################################## log2R boxplot 
######################################################## aneuploid / diploid cells
######################################################## experimental code
#####################################################################
#####################################################################

library(reshape2)
library(ggplot2)
library(dplyr)

# Extract logR matrix (columns 4 onward)
logR_mat <- CNA[, 4:ncol(CNA)]

# CopyKAT predictions
cell_status <- copykat_result$prediction
cell_status_vector <- cell_status$copykat.pred
names(cell_status_vector) <- cell_status$cell.names

# Normalize column names
colnames(logR_mat) <- gsub("\\.", "-", colnames(logR_mat))

# Keep only aneuploid and diploid
cell_status_vector <- cell_status_vector[cell_status_vector %in% c("aneuploid", "diploid")]

# Align predictions with CNA columns
cell_predictions <- cell_status_vector[colnames(logR_mat)]
valid_cells <- !is.na(cell_predictions)

logR_mat <- logR_mat[, valid_cells]
cell_predictions <- cell_predictions[valid_cells]

# Melt to long format
df_long <- melt(logR_mat, variable.name = "Cell", value.name = "logR")
df_long$Prediction <- cell_predictions[as.character(df_long$Cell)]

# Summary
summary_stats <- df_long %>%
  group_by(Prediction) %>%
  summarise(
    n = n(),
    min = min(logR, na.rm = TRUE),
    median = median(logR, na.rm = TRUE),
    max = max(logR, na.rm = TRUE)
  )

wilcox_test <- wilcox.test(logR ~ Prediction, data = df_long)

# Save stats
stats_file <- file.path(output_dir, "log2R_violin_stats_filtered.txt")
sink(stats_file)
cat("Summary Statistics (Filtered for Aneuploid/Diploid Only):\n\n")
print(summary_stats)
cat("\nWilcoxon Rank-Sum Test Result:\n")
print(wilcox_test)
sink()

# Violin plot with smaller font sizes
p_violin <- ggplot(df_long, aes(x = Prediction, y = logR, fill = Prediction)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5,
               fill = "white", color = "black") +
  scale_fill_manual(values = c("aneuploid" = "#D73027", "diploid" = "#4575B4")) +
  theme_bw() +
  labs(
    title = "SegMean by copykat class",
    x = "Prediction",
    y = "log2R"
  ) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "none"
  )

# Save plot
ggsave(
  filename = file.path(output_dir, "violin_log2R_by_prediction.png"),
  plot = p_violin,
  width = 4,
  height = 4,
  dpi = 300
)

#####################################################################
##################################################################### VIOLIN PLOT
######################################################## transform SegMean into Integer CN
######################################################## aneuploid / diploid cells
######################################################## experimental code 
#####################################################################
#####################################################################
##################################################################### COPY NUMBER
#####################################################################

# To transform SegMean into Abs Copy NUmber :
# R example
# seg.mean <- -0.5
# estimated_CN <- 2 * 2^seg.mean  # ~1.41   # a formula that different from usual WGS analysis
# round(estimated_CN)  # 1

# seg.mean  Estimated CN (rounded)  Interpretation
# 0.0 2.0 Normal (diploid)
# 0.58  ~3.0  Gain
# 1.0 4.0 High gain (e.g. amplification)
# -0.58 ~1.0  Single-copy loss
# -1.0  ~0.5  Strong loss / deletion

# Transform: logR to estimated integer CN
cn_estimated <- round(2 * 2^as.matrix(logR_mat))

# Convert back to data frame if needed
cn_estimated <- as.data.frame(cn_estimated)

# Preserve original column and row names
colnames(cn_estimated) <- colnames(logR_mat)
rownames(cn_estimated) <- rownames(logR_mat)

# Check structure
str(cn_estimated)

library(reshape2)
library(ggplot2)
library(dplyr)

# Melt to long format
df_long <- melt(cn_estimated, variable.name = "Cell", value.name = "copy_number")
df_long$Prediction <- cell_predictions[as.character(df_long$Cell)]

# Summary statistics by prediction group
summary_stats <- df_long %>%
  group_by(Prediction) %>%
  summarise(
    n = n(),
    min = min(copy_number, na.rm = TRUE),
    median = median(copy_number, na.rm = TRUE),
    max = max(copy_number, na.rm = TRUE)
  )

# Wilcoxon test
wilcox_test <- wilcox.test(copy_number ~ Prediction, data = df_long)

# Save summary + test
stats_file <- file.path(output_dir, "copy_number_violin_stats_filtered.txt")
sink(stats_file)
cat("Summary Statistics (Integer CN, Aneuploid/Diploid Only):\n\n")
print(summary_stats)
cat("\nWilcoxon Rank-Sum Test Result:\n")
print(wilcox_test)
sink()

# Optional: print to console
print(summary_stats)
print(wilcox_test)

# Violin plot using estimated CN values
p_violin <- ggplot(df_long, aes(x = Prediction, y = copy_number, fill = Prediction)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5,
               fill = "white", color = "black") +
  scale_fill_manual(values = c("aneuploid" = "#D73027", "diploid" = "#4575B4")) +
  theme_bw() +
  labs(
    title = "Estimated Copy Number by CopyKAT Class",
    x = "Prediction",
    y = "Estimated Copy Number"
  ) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "none"
  )

# Save the plot
ggsave(
  filename = file.path(output_dir, "violin_CN_by_prediction.png"),
  plot = p_violin,
  width = 4,
  height = 4,
  dpi = 300
)

#####################################################################
##################################################################### GISTIC2
######################################################## experimental code 
#####################################################################
#####################################################################

library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(ggplot2)
library(dplyr)

# https://github.com/bzhanglab/GISTIC2_example
# hg19.UCSC.add_miR.140312.refgene.mat
# hg38.UCSC.add_miR.160920.refgene.mat

# colnames(seg_df) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
# gistic2(
#  seg = "seg_file",  # Path to your SEG file (seq_file)
#  refgene = "/home/tanasa/hg38.gtf/hg38.UCSC.add_miR.160920.refgene.mat",
# )

#####################################################################
#####################################################################

# 🚨 Why Might Deletions Look Like Amplifications in CopyKAT?

# Cause Effect
# No proper normal reference : wrong log2 baseline
# Tumor-only input  Diploid cells misclassified
# Transcriptional noise Expression mismatch
# Smoothing artifacts Inversion of CNV signal
# Wrong genome version  Coordinate misalignment

#####################################################################
##################################################################### test gistic2

# https://github.com/bzhanglab/GISTIC2_example
# sample  chromosome      start   end     num.mark        log2

# devtools::install_github("WangLabCSU/blit")

# seg_df = seg
# colnames(seg_df) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

# seg_file

# colnames(seg_df) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
# blit::gistic2(
#  seg = "seg_file",  # Path to your SEG file (seq_file)
#  refgene = "/home/tanasa/hg38.gtf/hg38.UCSC.add_miR.160920.refgene.mat",
# )

# <Command: gistic2>