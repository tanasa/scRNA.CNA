
# Load required libraries

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
library(pheatmap)
library(proxy)
library(grid)

# Define marker sets (same as before)
marker_sets <- list(
  Astrocyte_like = c("GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A2", "GJA1"),
  Oligo_like     = c("OLIG1", "OLIG2", "SOX10", "PDGFRA", "PLP1", "CSPG4", "SOX6", "NKX2-2"),
  Mesenchymal    = c("CD44", "CHI3L1", "FN1", "SERPINE1", "VIM", "SPARC", "B2M", "CD74"),
  Progenitor     = c("SOX2", "NES", "DCX", "FABP7", "ASCL1", "HES6", "DLL3", "TUBB3"),
  Cycling        = c("MKI67", "TOP2A", "PCNA", "CDK1", "CENPF", "AURKB"),
  Glioma         = c("SOX2", "OLIG2", "EGFR", "PDGFRA", "S100A4", "NES", "CD44", "CHI3L1", "CDKN2A", "IDH1"),
  Oligodendrocytes = c("MBP", "PLP1", "MOG", "MAG", "SOX10", "OLIG1", "CNP", "GPR17", "TF", "CLDN11"),
  Pericytes      = c("PDGFRB", "RGS5", "ACTA2", "MYH11", "MCAM", "TAGLN", "NOTCH3"),
  Endothelial    = c("PECAM1", "VWF", "CD34", "CLDN5", "FLT1", "TIE1", "ENG"),
  B_Cells        = c("MS4A1", "CD19", "CD79A", "CD79B", "CD22", "CD24", "BLK"),
  T_Cells        = c("CD2", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "IL7R", "TRAC", "CCR7"),
  Myeloid        = c("CD14", "LYZ", "FCGR3A", "ITGAM", "CX3CR1", "C1QA", "AIF1", "TREM2", "P2RY12"),
  Macrophages    = c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"),
  NPC_like       = c("SOX2", "NES", "DLL3", "ASCL1", "HES6", "TUBB3", "FABP7", "DCX", "STMN2", "GPM6A")
)

# Get list of .rds files
rds_files <- list.files(pattern = "\\.rds$")

# Loop over each RDS file
for (rds_path in rds_files) {
  # Extract name without extension for labeling
  sample_name <- tools::file_path_sans_ext(basename(rds_path))
  message("Processing: ", sample_name)

  obj <- readRDS(rds_path)

  # Plot UMAP with cluster labels
  p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 1, label.size = 6)
  p$layers[[2]]$aes_params$fontface <- "bold"
  ggsave(filename = paste0(sample_name, "_UMAP_clusters.png"), plot = p, width = 8, height = 6, dpi = 300)

  # Create per-gene UMAP plots
  gene_umap_dir <- paste0(sample_name, "_marker_umaps")
  if (!dir.exists(gene_umap_dir)) dir.create(gene_umap_dir)

  for (cell_type in names(marker_sets)) {
    markers <- marker_sets[[cell_type]]
    available_genes <- markers[markers %in% rownames(obj)]

    for (gene in available_genes) {
      message("  Plotting gene: ", gene, " for ", cell_type)
      p <- FeaturePlot(obj, features = gene, reduction = "umap") +
        ggtitle(paste0(cell_type, ": ", gene)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))

      ggsave(file.path(gene_umap_dir, paste0(cell_type, "_", gene, ".png")),
             plot = p, width = 6, height = 5, dpi = 300)
    }
  }

  # Add module scores
  for (cell_type in names(marker_sets)) {
    gene_set <- marker_sets[[cell_type]]
    available_genes <- gene_set[gene_set %in% rownames(obj)]
    if (length(available_genes) > 0) {
      obj <- AddModuleScore(obj, features = list(available_genes), name = paste0(cell_type, "_Score"))
    }
  }

  # Plot module scores
  score_umap_dir <- paste0(sample_name, "_module_score")
  if (!dir.exists(score_umap_dir)) dir.create(score_umap_dir)

  module_score_cols <- grep("_Score1$", colnames(obj@meta.data), value = TRUE)
  for (score in module_score_cols) {
    message("  Plotting module score: ", score)
    p <- FeaturePlot(obj, features = score, reduction = "umap") +
      ggtitle(score) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8))

    ggsave(file.path(score_umap_dir, paste0("module_score_umaps_", score, ".png")),
           plot = p, width = 6, height = 5, dpi = 300)
  }

  # Optional: Save modified object
  saveRDS(obj, file = paste0(sample_name, "_scored.rds"))
}
