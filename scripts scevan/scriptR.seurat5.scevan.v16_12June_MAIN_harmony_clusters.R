######################################################################################
######################################################################################

# ----------------------------- #
#         Load Libraries        #
# ----------------------------- #

library(SCEVAN)
library(Seurat)
library(dplyr)
library(SeuratData)
library(SeuratObject)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(tibble)
library(dplyr)
library(readr)
library(forcats)
# library(GenomicPlot)
# library(karyoploteR)
library(GenomicRanges)
library(ggplot2)
library(circlize)
library(proxy)
library(parallelDist)
library(RColorBrewer)
library(grid)

# Suppress ComplexHeatmap messages
library(ComplexHeatmap) 
ht_opt$message = FALSE

######################################################################################
######################################################################################

options(future.globals.maxSize = 1e9)

# For Unix/Linux systems, use multicore instead of SOCK
library(BiocParallel)
param <- MulticoreParam(workers = 8)  # Reduce workers
register(param)

# Sys.setenv("BIOCPARALLEL_WORKERS" = "1")
# library(BiocParallel)
# register(SerialParam())

######################################################################################
######################################################################################

# ----------------------------- #
#         Load RDS              #
# ----------------------------- #

# Read file and extract clean filename
filename <- "LGG-03.rds"
obj <- readRDS(filename)

# Extract filename without extension for pasting
file_base <- gsub("\\.rds$", "", basename(filename))
# Result: "copykat_obj_sub"

# ----------------------------- #
#     Setup Output Directory    #
# ----------------------------- #

# method_name <- "harmony"
method_name = ""

# Create full output path including "output/" prefix that SCEVAN expects
# scevan_base_dir <- file.path("output", paste0(method_name, "scevan_output"))
# scevan_dir <- file.path(scevan_base_dir, method_name)
# scevan_dir = "output"

# Now you can paste with other words
scevan_dir <- paste(file_base, "output", sep = ".")

dir.create(scevan_dir, recursive = TRUE, showWarnings = FALSE)

# Define sample output prefix
# sample_prefix <- file.path(scevan_dir, "scevan")
# dir.create(dirname(sample_prefix), recursive = TRUE, showWarnings = FALSE)

cat("Running SCEVAN for method:", method_name, "\n")
cat("Output directory:", scevan_dir, "\n\n")

# setwd("./output/")

######################################################################################
######################################################################################

# ----------------------------- #
#     Prepare Seurat Object     #
# ----------------------------- #

obj_scevan <- obj  # Clone original Seurat object

# str(obj@assays)
# str(obj@assays$RNA)
# str(obj@assays$RNA@layers)
# str(obj@assays$RNA@layers$data)
# str(obj@assays$RNA@layers$counts)
# str(obj@meta.data)

cat("Original RNA layers:\n")
print(Layers(obj_scevan[["RNA"]]))

obj_scevan <- JoinLayers(obj_scevan)

cat("After JoinLayers RNA layers:\n")
print(Layers(obj_scevan[["RNA"]]))

cat("\nStructure of joined object:\n")
str(obj_scevan)

cat("\nDimensions of layers in original object:\n")
for (layer_name in Layers(obj[["RNA"]])) {
  dims <- dim(obj[["RNA"]]@layers[[layer_name]])
  cat(sprintf("  %-25s : %d x %d\n", layer_name, dims[1], dims[2]))
}

cat("\nDimensions of layers in joined object:\n")
for (layer_name in Layers(obj_scevan[["RNA"]])) {
  dims <- dim(obj_scevan[["RNA"]]@layers[[layer_name]])
  cat(sprintf("  %-25s : %d x %d\n", layer_name, dims[1], dims[2]))
}

######################################################################################
###################################################################################### UMAP visualization

# Run PCA (if not done already)
obj_scevan <- RunPCA(obj_scevan, assay = "RNA", features = VariableFeatures(obj_scevan))

# Run UMAP
obj_scevan <- RunUMAP(obj_scevan, dims = 1:20, reduction = "pca", assay = "RNA")

# Generate the plot with larger points
umap_plot <- DimPlot(obj_scevan, reduction = "umap", 
                     group.by = "seurat_clusters", 
                     pt.size = 1.5) + 
                     ggtitle(file_base)

# Display plot
print(umap_plot)

# Save plot to file
output_file <- file.path(scevan_dir, paste0(file_base,"_seurat_object.pdf"))
ggsave(output_file, plot = umap_plot, width = 8, height = 6)

######################################################################################
######################################################################################
# ----------------------------- #
#         Run SCEVAN            #
# ----------------------------- #
######################################################################################
######################################################################################

cat("Running SCEVAN pipeline...\n")
count_mtx <- as.matrix(LayerData(obj_scevan[["RNA"]], layer = "counts"))

scevan_results <- pipelineCNA(
  count_mtx = count_mtx,
  # sample = file.path(scevan_dir, "scevan"),  # Prefix for output files
  sample = "",        # Sample name to save results (optional)
  par_cores = 10,     # Number of cores to run the pipeline (optional - default 20)
  norm_cell = NULL,
  SUBCLONES = TRUE,
  beta_vega = 0.5,
  ClonalCN = TRUE,
  plotTree = TRUE,
# AdditionalGeneSets = NULL,
# SCEVANsignatures = TRUE,
  organism = "human",
  ngenes_chr = 5,
  perc_genes = 10,
  FIXED_NORMAL_CELLS = FALSE
)

# Save RDS file
saveRDS(scevan_results, 
        file = file.path(scevan_dir, "scevan_results.rds"))

# Save CSV file
write.csv(scevan_results, 
          file = file.path(scevan_dir, "scevan_results.csv"), 
          row.names = FALSE)

# Or if it's a list of data frames
# lapply(names(scevan_results), function(name) {
#       write.csv(scevan_results[[name]], file = paste0("scevan_results_", name, ".csv"), row.names = FALSE)
# })

# pipelineCNA Executes the entire SCEVAN pipeline that classifies
#     tumour and normal cells from the raw count matrix, infer the
#     clonal profile of cancer cells and looks for possible sub-clones
#     in the tumour cell matrix automatically analysing the specific and
#     shared alterations of each subclone and a differential analysis of
#     pathways and genes expressed in each subclone.

# Usage:
#
#     pipelineCNA(
#       count_mtx,
#       sample = "",
#       par_cores = 20,
#       norm_cell = NULL,
#       SUBCLONES = TRUE,
#       beta_vega = 0.5,
#       ClonalCN = TRUE,
#       plotTree = TRUE,
#       AdditionalGeneSets = NULL,
#       SCEVANsignatures = TRUE,
#       organism = "human",
# 
#       ngenes_chr = 5,
#       perc_genes = 10,
#       FIXED_NORMAL_CELLS = FALSE
#     )

# SUBCLONES: Boolean value TRUE if you are interested in analysing the
#          clonal structure and FALSE if you are only interested in the
#          classification of malignant and non-malignant cells (optional
#          - default TRUE)
#
# beta_vega: Specifies beta parameter for segmentation, higher beta for
#          more coarse-grained segmentation. (optional - default 0.5)
#
# ClonalCN: Get clonal CN profile inference from all tumour cells
#          (optional)
#
# plotTree: Plot Phylogenetic tree (optional - default FALSE)
# 
# AdditionalGeneSets: list of additional signatures of normal cell types
#          (optional)
#
# SCEVANsignatures: FALSE if you only want to use only the signatures
#          specified in AdditionalGeneSets (default TRUE)
#
#          FIXED_NORMAL_CELLS: TRUE if norm_cell vector to be used as fixed
#          reference, if you are only interested in clonal structure and
#          not normal/tumor classification (default FALSE)
#
# Examples:
#
#     res_pip <- pipelineCNA(count_mtx)
#
#
# head(scevan_results)
#                                       class confidentNormal subclone
# LGG-04-1_LGG-04-1_AAACCCAAGTCCCAGC-1 normal            <NA>       NA
# LGG-04-1_LGG-04-1_AAACCCAAGTGCCGAA-1 normal            <NA>       NA
# LGG-04-1_LGG-04-1_AAACCCACAAAGTGTA-1  tumor            <NA>        1
# LGG-04-1_LGG-04-1_AAACCCACACTGGAAG-1  tumor            <NA>        2
# LGG-04-1_LGG-04-1_AAACCCATCGACCATA-1  tumor            <NA>        3

# careful with : step 9: Segmentation, where the function tries to construct a data
# This issue is related to instability in parallel::mclapply() when used with many 
# cores or on certain systems (especially non-Linux or shared environments).
# Try setting par_cores = 1 to force sequential execution

# head(scevan_results)
#                                         class confidentNormal subclone
#
#                                         class confidentNormal subclone
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1 filtered            <NA>       NA
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1   normal            <NA>       NA

# table(scevan_results$subclone)
# 
#  1  2  3  4  5 
# 62 65 90 94 79 

# table(scevan_results$confidentNormal)
# yes 
# 30

######################################################################################
######################################################################################
######################################################################################  RESULTS

# the resulting files are the following :

# Clonal_CN.seg Overall tumor segmentation (GISTIC-style)
# CNAmtx.RData  Matrix of CNA values across all cells
# CNAmtxSubclones.RData CNA profiles per subclone
# count_mtx_annot.RData Annotated expression data
# onlytumorvega_output  Interactive visualization for tumor cells
# PlotOncoHeat.RData  CNA heatmap object
# subclone[1–5]_CN.seg  CNA profiles per subclone
# SubcloneDiffAnalysis.RData  Stat analysis of CNA differences between clones
# vega_output Interactive CNA plot (full dataset)

setwd("./output/")
getwd()

######################################################################################
######################################################################################
######################################################################################
######################################################################################

# ✅ Find the number of .seg files in the output directory
# ✅ Read them into R as a named list of data frames
# ✅ Build a data structure (seg_files_data) where each entry corresponds to the subclone/clonal .seg file


# List all .seg files in the directory
all_seg_files <- list.files("./", pattern = "\\.seg$", full.names = TRUE)

# Print how many .seg files were found
cat("Number of .seg files:", length(all_seg_files), "\n")

# List all .seg files in the directory
# Define expected names for subclone files
# seg_files <- list(
#  Clonal     = "_Clonal_CN.seg",
#  Subclone1  = "_subclone1_CN.seg",
#  Subclone2  = "_subclone2_CN.seg",
#  Subclone3  = "_subclone3_CN.seg",
#  Subclone4  = "_subclone4_CN.seg",
#  Subclone5  = "_subclone5_CN.seg"
# )

# Now:
# ➜ Each data frame is accessible as 'clonal', 'subclone1', 'subclone2', ...
# ➜ seg_files_data is a list holding all of them for easy access!

# Initialize a list to hold data frames
seg_files_data <- list()

# Loop through each .seg file
for (file_path in all_seg_files) {
  # Get the filename (without path)
  file_name <- basename(file_path)
  
  # Determine the subclone name (e.g., subclone1, subclone2, etc.) or "clonal"
  if (grepl("Clonal", file_name, ignore.case = TRUE)) {
    var_name <- "clonal"
  } else if (grepl("subclone[0-9]+", file_name, ignore.case = TRUE)) {
    subclone_num <- regmatches(file_name, regexpr("subclone[0-9]+", file_name, ignore.case = TRUE))
    var_name <- tolower(subclone_num)  # e.g., subclone1
  } else {
    next  # skip files that don't match these patterns
  }
  
  # Read the file into a data frame
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Store the data frame in the list
  seg_files_data[[var_name]] <- df
  
  cat("Loaded:", var_name, "from", file_name, "\n")
}

# Now seg_files_data is a named list of data frames
# Check what's inside

str(seg_files_data)

# Find how many data frames (clones) are stored in seg_files_data
num_dataframes <- length(seg_files_data)
cat("Number of main clones and subclones in seg_files_data:", num_dataframes, "\n")

######################################################################################
###################################################################################### customized visualization
###################################################################################### display the seg mean
###################################################################################### 

# Chr Pos End CN  segm.mean
# 1 1 944204  11806920  0 -0.206613
# 2 1 11806096  30757820  0 -0.44116
# 3 1 30869467  35557950  1 -0.147513
# 4 1 35808172  51273455  0 -0.371267

# Define segmentation files
# seg_files <- list(
#  Clonal     = "_Clonal_CN.seg",
#  Subclone1  = "_subclone1_CN.seg",
#  Subclone2  = "_subclone2_CN.seg",
#  Subclone3  = "_subclone3_CN.seg",
#  Subclone4  = "_subclone4_CN.seg",
#  Subclone5  = "_subclone5_CN.seg"
# )

# Create a list of data frames adding clone name
seg_list <- lapply(names(seg_files_data), function(clone_name) {
  df <- seg_files_data[[clone_name]]
  colnames(df) <- c("chr", "start", "end", "cn", "segm.mean")
  df$clone <- clone_name
  return(df)
})

# Combine into one big data frame
seg_all <- bind_rows(seg_list)

# Clean chromosome names
seg_all$chr <- gsub("^chr", "", seg_all$chr)
seg_all$chr <- factor(seg_all$chr, levels = c(1:22, "X", "Y"))

# Consistent ordering of clones
seg_all$clone <- factor(seg_all$clone, levels = rev(names(seg_files_data)))

# Plot each chromosome
plot_list <- lapply(levels(seg_all$chr), function(chrname) {
  df_chr <- filter(seg_all, chr == chrname)
  
  ggplot(df_chr, aes(xmin = start, xmax = end,
                     ymin = as.numeric(clone) - 0.4, ymax = as.numeric(clone) + 0.4,
                     fill = segm.mean)) +
    geom_rect(color = "white") +
    scale_y_continuous(breaks = seq_along(levels(seg_all$clone)),
                       labels = rev(levels(seg_all$clone))) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "segm.mean", limits = c(-0.5, 0.5)) +
    labs(title = paste("Chromosome", chrname), x = "Genome Position (bp)", y = "Clone") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
})

# Save to PDF
pdf("genome_wide_segm_mean_heatmaps_all_chromosomes.pdf", width = 16, height = 40)
wrap_plots(plot_list, ncol = 2)
dev.off()

# Combine and save to PNG
combined_plot <- wrap_plots(plot_list, ncol = 2)

ggsave("genome_wide_segm_mean_heatmaps_all_chromosomes.png",
       combined_plot,
       width = 16, height = 40, dpi = 300, units = "in")

######################################################################################
###################################################################################### customized visualization
###################################################################################### display the CN calls
###################################################################################### 

# seg_files <- list(
#  Clonal     = "_Clonal_CN.seg",
#  Subclone1  = "_subclone1_CN.seg",
#  Subclone2  = "_subclone2_CN.seg",
#  Subclone3  = "_subclone3_CN.seg",
#  Subclone4  = "_subclone4_CN.seg",
#  Subclone5  = "_subclone5_CN.seg"
# )

# Create a list of data frames adding clone name
seg_list <- lapply(names(seg_files_data), function(clone_name) {
  df <- seg_files_data[[clone_name]]
  colnames(df) <- c("chr", "start", "end", "cn", "segm.mean")
  df$clone <- clone_name
  return(df)
})

# Combine into one big data frame
seg_all <- bind_rows(seg_list)

# Clean chromosome names
seg_all$chr <- gsub("^chr", "", seg_all$chr)
seg_all$chr <- factor(seg_all$chr, levels = c(1:22, "X", "Y"))

# Consistent ordering of clones
seg_all$clone <- factor(seg_all$clone, levels = rev(names(seg_files_data)))

# Plot each chromosome with CN numerical values
plot_list <- lapply(levels(seg_all$chr), function(chrname) {
  df_chr <- filter(seg_all, chr == chrname)
  
  ggplot(df_chr, aes(xmin = start, xmax = end,
                     ymin = as.numeric(clone) - 0.4, ymax = as.numeric(clone) + 0.4,
                     fill = cn)) +
    geom_rect(color = "white") +
    geom_text(aes(x = (start + end)/2, 
                   y = as.numeric(clone), 
                   label = round(cn, 1)), 
               size = 2.5, color = "black") +  # Adjust size/color as needed
    scale_y_continuous(breaks = seq_along(levels(seg_all$clone)),
                       labels = rev(levels(seg_all$clone))) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2,
                         name = "CN", limits = c(0, 4)) +
    labs(title = paste("Chromosome", chrname), x = "Genomic Position (bp)", y = "Clone") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
})

# Save to file:
pdf("genome_wide_copy_numbers_heatmaps_all_chromosomes.pdf", width = 16, height = 40)
wrap_plots(plot_list, ncol = 2)
dev.off()

# Combine and save to PNG
combined_plot <- wrap_plots(plot_list, ncol = 2)

ggsave("genome_wide_copy_numbers_heatmaps_all_chromosomes.png",
       combined_plot,
       width = 16, height = 40, dpi = 300, units = "in")


######################################################################################
######################################################################################
######################################################################################

# seg_files <- list(
#  Clonal     = "_Clonal_CN.seg",
#  Subclone1  = "_subclone1_CN.seg",
#  Subclone2  = "_subclone2_CN.seg",
#  Subclone3  = "_subclone3_CN.seg",
#  Subclone4  = "_subclone4_CN.seg",
#  Subclone5  = "_subclone5_CN.seg"
# )

# Create a list of data frames adding clone name
seg_list <- lapply(names(seg_files_data), function(clone_name) {
  df <- seg_files_data[[clone_name]]
  colnames(df) <- c("chr", "start", "end", "cn", "segm.mean")
  df$clone <- clone_name
  return(df)
})

# Combine seg_list into a single data frame
seg_all <- bind_rows(seg_list)

# Format chromosome and clone factors
seg_all$chr <- gsub("^chr", "", seg_all$chr)
chrom_levels <- as.character(c(1:22, "X", "Y"))
seg_all <- seg_all %>% filter(chr %in% chrom_levels)
seg_all$chr <- factor(seg_all$chr, levels = chrom_levels)
seg_all$clone <- factor(seg_all$clone, levels = rev(names(seg_files_data)))  # use actual loaded data

# Build the plot showing segm.mean
p <- ggplot(seg_all, aes(xmin = start, xmax = end,
                         ymin = as.numeric(clone) - 0.4,
                         ymax = as.numeric(clone) + 0.4,
                         fill = segm.mean)) +  # Fill based on segm.mean
  geom_rect(color = NA) +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
    values = scales::rescale(c(-0.5, -0.2, 0, 0.2, 0.5)),
    limits = c(-0.5, 0.5),
    name = "Segm. Mean"
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(seg_all$clone)),
    labels = levels(seg_all$clone),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Segment Mean per Clone per Chromosome",
    x = NULL,
    y = "Clone"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

# Save the plot with a shorter vertical height
ggsave("genome_wide_segm_means_heatmaps_per_chromosome.pdf",
       plot = p, width = 20, height = 2, dpi = 300)
ggsave("genome_wide_segm_means_heatmaps_per_chromosome.png",
       plot = p, width = 20, height = 2, dpi = 300, bg = "white")

######################################################################################
######################################################################################

# seg_files <- list(
#  Clonal     = "_Clonal_CN.seg",
#  Subclone1  = "_subclone1_CN.seg",
#  Subclone2  = "_subclone2_CN.seg",
#  Subclone3  = "_subclone3_CN.seg",
#  Subclone4  = "_subclone4_CN.seg",
#  Subclone5  = "_subclone5_CN.seg"
# )

# Create a list of data frames adding clone name
seg_list <- lapply(names(seg_files_data), function(clone_name) {
  df <- seg_files_data[[clone_name]]
  colnames(df) <- c("chr", "start", "end", "cn", "segm.mean")
  df$clone <- clone_name
  return(df)
})

# Combine seg_list into a single data frame
seg_all <- bind_rows(seg_list)

# Format chromosome and clone factors
seg_all$chr <- gsub("^chr", "", seg_all$chr)
chrom_levels <- as.character(c(1:22, "X", "Y"))
seg_all <- seg_all %>% filter(chr %in% chrom_levels)
seg_all$chr <- factor(seg_all$chr, levels = chrom_levels)
seg_all$clone <- factor(seg_all$clone, levels = rev(names(seg_files_data)))  # use actual loaded data

# Build the plot
p <- ggplot(seg_all, aes(xmin = start, xmax = end,
                         ymin = as.numeric(clone) - 0.4,
                         ymax = as.numeric(clone) + 0.4,
                         fill = cn)) +  # Change to segm.mean if needed!
  geom_rect(color = NA) +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
    values = scales::rescale(c(0, 1, 2, 3, 4)),  # Adjust for CN range
    limits = c(min(seg_all$cn, na.rm = TRUE), max(seg_all$cn, na.rm = TRUE)),
    name = "Copy Number"
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(seg_all$clone)),
    labels = levels(seg_all$clone),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Copy Number",
    x = NULL,
    y = "Clone"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

# Save the plot with a shorter vertical height
ggsave("genome_wide_copy_numbers_heatmaps_per_chromosome.pdf",
       plot = p, width = 20, height = 2, dpi = 300) 
ggsave("genome_wide_copy_numbers_heatmaps_per_chromosome.png",
       plot = p, width = 20, height = 2, dpi = 300, bg = "white")

######################################################################################
######################################################################################
######################################################################################
######################################################################################

# setwd("./output/")
# reading the MATRIX files provided by SCEVAN :

# _CNAmtxSubclones.RData
# _count_mtx_annot.RData

######################################################################################
######################################################################################
######################################################################################
######################################################################################

# We use new.env() with load(..., envir = env4) to isolate the objects that are loaded from the .RData file. 

# env4 <- new.env()
# load("file.RData", envir = env4)

# Access objects selectively
# get("oncoHeat", envir = env4)

######################################################################################
######################################################################################
###################################################################################### working with : 
###################################################################################### CNA matrices
# -----------------------------
# 1. Load _CNAmtx.RData
# -----------------------------
###################################################################################### env1
######################################################################################
###################################################################################### the object is : 
###################################################################################### "CNA_mtx_relat"

# cat("\n--- Loading: _CNAmtx.RData ---\n")
# env1 <- new.env()
# load("_CNAmtx.RData", envir = env1)
# print(ls(env1))

# for (obj_name in ls(env1)) {
#  cat("\nObject:", obj_name, "\n")
#  obj <- get(obj_name, envir = env1)
#  print(class(obj))

#  if (is.matrix(obj) || is.data.frame(obj)) {
#    print(dim(obj))
#    cat("Column names:\n")
#    print(head(colnames(obj), 3))
#    cat("Preview [3x3]:\n")
#    print(obj[1:min(3, nrow(obj)), 1:min(3, ncol(obj))])
#  } else if (is.list(obj)) {
#    print(str(obj, max.level = 1))
#  } else {
#    print(obj)
#  }
# }

# -- Loading: _CNAmtx.RData ---
# [1] "CNA_mtx_relat"

# Object: CNA_mtx_relat 
# [1] "matrix" "array" 
# [1] 7166  952
# Column names:
# [1] "LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1"
# [2] "LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1"
# [3] "LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1"
# Preview [3x3]:
#     LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1 LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1
# [1,]                           -0.1706296                          -0.02180575
# [2,]                           -0.1706296                          -0.02180575
# [3,]                           -0.1706296                          -0.02180575
#     LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1
# [1,]                            0.1550459
# [2,]                            0.1550459
# [3,]                            0.1550459

load("_CNAmtx.RData")
str(CNA_mtx_relat)

######################################################################################
######################################################################################
###################################################################################### working with : 
###################################################################################### CNA matrices
# -----------------------------
# 2. Load _CNAmtxSubclones.RData
# -----------------------------
###################################################################################### the object is : 
###################################################################################### results.com
###################################################################################### 
###################################################################################### env2

# cat("\n--- Loading: _CNAmtxSubclones.RData ---\n")
# env2 <- new.env()
# load("_CNAmtxSubclones.RData", envir = env2)
# print(ls(env2))

# for (obj_name in ls(env2)) {
#  cat("\nObject:", obj_name, "\n")
#  obj <- get(obj_name, envir = env2)
#  print(class(obj))

#  if (is.matrix(obj) || is.data.frame(obj)) {
#    print(dim(obj))
#    cat("Column names:\n")
#    print(head(colnames(obj), 3))
#    cat("Preview [3x3]:\n")
#    print(obj[1:min(3, nrow(obj)), 1:min(3, ncol(obj))])
#  } else if (is.list(obj)) {
#    print(str(obj, max.level = 1))
#  } else {
#    print(obj)
#  }
# }

# Object: results.com 
# [1] "matrix" "array" 
# [1] 7166  390
# Column names:
# [1] "LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1"
# [2] "LGG-04-1_LGG-04-1_AATAGAGAGACCTCAT-1"
# [3] "LGG-04-1_LGG-04-1_ACGGGTCCACCACATA-1"
# Preview [3x3]:
#     LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1 LGG-04-1_LGG-04-1_AATAGAGAGACCTCAT-1
# [1,]                           -0.2358805                          -0.09943377
# [2,]                           -0.2358805                          -0.09943377
# [3,]                           -0.2358805                          -0.09943377
#      LGG-04-1_LGG-04-1_ACGGGTCCACCACATA-1
# [1,]                          -0.09724475
# [2,]                          -0.09724475
# [3,]                          -0.09724475

load("_CNAmtxSubclones.RData")
CNA_mtx_subclones = results.com 
str(CNA_mtx_subclones)

######################################################################################
######################################################################################
###################################################################################### the object is : 
# -----------------------------
# 3. Load _count_mtx_annot.RData
# -----------------------------
###################################################################################### _count_mtx_annot
######################################################################################
###################################################################################### env3

# cat("\n--- Loading: _count_mtx_annot.RData ---\n")
# env3 <- new.env()
# load("_count_mtx_annot.RData", envir = env3)
# print(ls(env3))

# for (obj_name in ls(env3)) {
#  cat("\nObject:", obj_name, "\n")
#  obj <- get(obj_name, envir = env3)
#  print(class(obj))

#  if (is.matrix(obj) || is.data.frame(obj)) {
#    print(dim(obj))
#    cat("Column names:\n")
#    print(head(colnames(obj), 5))
#    cat("Preview [3x3]:\n")
#    print(obj[1:min(5, nrow(obj)), 1:min(5, ncol(obj))])
#  } else if (is.list(obj)) {
#    print(str(obj, max.level = 1))
#  } else {
#    print(obj)
#  }
# }

# Object: count_mtx_annot 
# [1] "data.frame"
# [1] 7166    5
# Column names:
# [1] "seqnames"  "start"     "end"       "gene_id"   "gene_name"
# Preview [3x3]:
#                seqnames   start     end         gene_id gene_name
# ENSG00000188976        1  944204  959309 ENSG00000188976     NOC2L
# ENSG00000188290        1  998962 1000172 ENSG00000188290      HES4
# ENSG00000187608        1 1001138 1014541 ENSG00000187608     ISG15
# ENSG00000078808        1 1216908 1232031 ENSG00000078808      SDF4
# ENSG00000176022        1 1232265 1235041 ENSG00000176022   B3GALT6

load("_count_mtx_annot.RData")
str(count_mtx_annot)

###################################################################################### the objects are :
###################################################################################### oncoHeat
###################################################################################### plotHeatmapOncoHeat
# -----------------------------
# 4. Load PlotOncoHeat.RData
# -----------------------------
######################################################################################
######################################################################################
###################################################################################### env4

# cat("\n--- Loading: PlotOncoHeat.RData ---\n")
# env4 <- new.env()
# load("PlotOncoHeat.RData", envir = env4)
# print(ls(env4))

# for (obj_name in ls(env4)) {
#  cat("\nObject:", obj_name, "\n")
#  obj <- get(obj_name, envir = env4)
#  print(class(obj))
#
#  if (is.matrix(obj) || is.data.frame(obj)) {
#    print(dim(obj))
#    cat("Column names:\n")
#    print(head(colnames(obj), 3))
#    cat("Preview [3x3]:\n")
#    print(obj[1:min(3, nrow(obj)), 1:min(3, ncol(obj))])
#  } else if (is.list(obj)) {
#    print(str(obj, max.level = 1))
#  } else {
#    print(obj)
#  }
# }

# --- Loading: PlotOncoHeat.RData ---
# [1] "oncoHeat"            "plotHeatmapOncoHeat"
# Object: oncoHeat 
# [1] "data.frame"
# [1]   5 109
# Column names:
# [1] "1 (p36.13-p33) "    "1 (p36.22-p35.2) "  "1 (p36.13-p36.11) "
# Preview [3x3]:
#          1 (p36.13-p33)  1 (p36.22-p35.2)  1 (p36.13-p36.11) 
# Subclone1               0                 0                  0
# Subclone2              -2                 0                  0
# Subclone3               0                 0                  0

# for (obj_name in ls(env4)) {
#  cat("\nObject:", obj_name, "\n")
#  obj <- get(obj_name, envir = env4)
#  print(class(obj)) 
# }

# Object: oncoHeat 
# [1] "data.frame"
# Object: plotHeatmapOncoHeat 
# [1] "pheatmap"

# oncoHeat <- get("oncoHeat", envir = env4)

# dim(oncoHeat)
# head(colnames(oncoHeat) ,2)
# head(oncoHeat, 2)

# [1]   5 109
# [1] "1 (p36.13-p33) "   "1 (p36.22-p35.2) "
#           1 (p36.13-p33)  1 (p36.22-p35.2)  1 (p36.13-p36.11) 
# Subclone1               0                 0                  0
# Subclone2              -2                 0                  0

# load("PlotOncoHeat.RData", envir = env4)

load("PlotOncoHeat.RData")
str(oncoHeat)

source("../script.plot.cytoband.R")

######################################################################################
###################################################################################### working with
######################################################################################
######################################################################################
######################################################################################
###################################################################################### VEGA

# Vega : 
#
# 🧬 Key Columns:
# Column  Description
# Chr, Start, End Genomic coordinates of each segment
# Size  Size of the segment (in bp)
# Mean  Average signal (likely log2 ratio or segm.mean)
# L.pv, G.pv, LOH.pv  P-values for loss, gain, and LOH
# X..L, X.G, X.LOH  % of the segment showing loss/gain/LOH
# Probe Size  Number of probes in the segment
# Loss Mean, Gain Mean, LOH Mean  Segment means conditioned on type
# Focal-score Loss, Focal-score Gain, Focal-score LOH Scores for focal events

# setwd("../")

######################################################################################
######################################################################################
###################################################################################### VEGA

# file1 <- "onlytumorvega_output"
# Read the file
# cnv_data1 <- read.delim(file1, stringsAsFactors = FALSE)

# Show basic info
# cat("Number of rows:", nrow(cnv_data1), "\n")
# cat("Number of columns:", ncol(cnv_data1), "\n")
# cat("Column names:\n")
# print(colnames(cnv_data1))

#  Chr     Start       End     Size      Mean  L.pv G.pv LOH.pv  X..L  X.G X.LOH
# 1   1    959309  19251552 18292244  0.066466 1.000    1      1  0.1% 0.2%    0%
# 2   1  19260128  28582983  9322856 -0.098005 0.000    1      1 18.9% 0.1%    0%
# 3   1  28643085  36464437  7821353  0.048705 0.991    1      1  0.3% 1.4%    0%
# 4   1  36483278  47232220 10748943 -0.020958 0.000    1      1  7.5% 0.2%    0%
# 5   1  47378839 109275750 61896912  0.095023 1.000    1      1    0% 4.2%    0%
# 6   1 109283186 114153919  4870734 -0.075478 0.000    1      1 16.2% 0.2%    0%
#  Probe.Size Loss.Mean Gain.Mean LOH.Mean Focal.score.Loss Focal.score.Gain
# 1         81 -0.224725  0.205858        0     2.913102e-06     5.339601e-06
# 2         53 -0.296726  0.202558        0     1.058562e-03     4.012942e-06
# 3         43 -0.216486  0.217654        0     1.586389e-05     6.911780e-05
# 4         56 -0.261920  0.208660        0     3.488213e-04     7.828476e-06
# 5        111  0.000000  0.223039        0     0.000000e+00     8.442729e-05
# 6         28 -0.313782  0.218660        0     1.812819e-03     1.640731e-05
#   Focal.score.LOH
# 1               0
# 2               0
# 3               0
# 4               0
# 5               0
# 6               0

# file2 <- "vega_output"
# Read the file
# cnv_data2 <- read.delim(file2, stringsAsFactors = FALSE)

# Show basic info
# cat("Number of rows:", nrow(cnv_data2), "\n")
# cat("Number of columns:", ncol(cnv_data2), "\n")
# cat("Column names:\n")
# print(colnames(cnv_data2))

# head(cnv_data2, 3)
#   Chr    Start      End     Size      Mean  L.pv G.pv LOH.pv  X..L  X.G X.LOH
# 1   1   959309 19251552 18292244  0.066466 0.999    1      1  0.1% 0.2%    0%
# 2   1 19260128 28582983  9322856 -0.098005 0.000    1      1 18.9% 0.1%    0%
# 3   1 28643085 36464437  7821353  0.048705 0.991    1      1  0.3% 1.4%    0%
#  Probe.Size Loss.Mean Gain.Mean LOH.Mean Focal.score.Loss Focal.score.Gain
# 1         81 -0.224725  0.205858        0     2.913102e-06     5.339601e-06
# 2         53 -0.296726  0.202558        0     1.058562e-03     4.012942e-06
# 3         43 -0.216486  0.217654        0     1.586389e-05     6.911780e-05
#  Focal.score.LOH
# 1               0
# 2               0
# 3               0

file1 <- "onlytumorvega_output"
file2 <- "vega_output"

# Read data frames
cnv_data1 <- read.delim(file1, stringsAsFactors = FALSE)
cnv_data2 <- read.delim(file2, stringsAsFactors = FALSE)

# Show basic info
cat("File 1 - Rows:", nrow(cnv_data1), "Columns:", ncol(cnv_data1), "\n")
cat("File 2 - Rows:", nrow(cnv_data2), "Columns:", ncol(cnv_data2), "\n")

cat("File 1 Column Names:\n")
print(colnames(cnv_data1))

cat("File 2 Column Names:\n")
print(colnames(cnv_data2))

# Head of file 1
head(cnv_data1, 3)
# Head of file 2
head(cnv_data2, 3)

# Dim of file 1
dim(cnv_data1)
# Dim of file 2
dim(cnv_data2)

# BOXPLOTS : Columns to plot

cols_to_plot <- c("Mean", "L.pv", "G.pv", "LOH.pv", "X..L", "X.G", "X.LOH",
                   "Probe.Size", "Loss.Mean", "Gain.Mean", "LOH.Mean",
                   "Focal.score.Loss", "Focal.score.Gain")

# Ensure these columns exist
common_cols <- intersect(cols_to_plot, names(cnv_data1))

# Prepare data in long format
df1_long <- data.frame(Value = as.vector(as.matrix(cnv_data1[, common_cols])),
                       Variable = rep(common_cols, each = nrow(cnv_data1)),
                       Source = "onlytumor")  # relabeled

df2_long <- data.frame(Value = as.vector(as.matrix(cnv_data2[, common_cols])),
                       Variable = rep(common_cols, each = nrow(cnv_data2)),
                       Source = "allcells")  # relabeled

# Combine long data frames
combined_long <- rbind(df1_long, df2_long)
combined_long$Value <- as.numeric(combined_long$Value)

# Create ggplot boxplots with 3 plots per row and improved aesthetics (NO jitter)
p <- ggplot(combined_long, aes(x = Source, y = Value, fill = Source)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "grey30", width = 0.4) +  # narrow boxplots, no outliers
  facet_wrap(~ Variable, scales = "free", nrow = ceiling(length(common_cols) / 3)) +  # 3 per row
  theme_bw(base_size = 14) +  # Base font size
  labs(y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("onlytumor" = "#66c2a5", "allcells" = "#fc8d62"))  # Custom fill colors

# Save the plot
ggsave("vega_param_comparison_onlytumor_vs_vegaoutput.pdf", plot = p, width = 6, height = 14)
ggsave("vega_param_comparison_onlytumor_vs_vegaoutput.png", plot = p, width = 6, height = 14, dpi = 300)

######################################################################################
######################################################################################
######################################################################################
###################################################################################### META DATA
######################################################################################
######################################################################################
######################################################################################

# ----------------------------- #
#     Add Metadata & Save       #
# ----------------------------- #

scevan_results2 <- scevan_results[colnames(obj_scevan), , drop = FALSE]
obj_scevan2 <- AddMetaData(obj_scevan, metadata = scevan_results2)

cat("\nSCEVAN metadata structure:\n")
head(scevan_results2, 2)

# head(colnames(obj_scevan2@meta.data),20)
# [1] "orig.ident"       "nCount_RNA"       "nFeature_RNA"     "sample"          
# [5] "percent.mt"       "RNA_snn_res.0.5"  "seurat_clusters"  "harmony_clusters"
# [9] "class"            "confidentNormal"  "subclone"

# Save the SCEVAN results to a CSV file
write.csv(scevan_results2, file = "scevan_results.csv", row.names = TRUE)
# write.table(scevan_results2, file = file.path(scevan_dir, "scevan_results2.txt"), row.names = TRUE)

# If you want to save the complete results (including all aspects of the SCEVAN analysis)
saveRDS(scevan_results, file = "scevan_results2_complete.rds")
# saveRDS(scevan_results, file = file.path(scevan_dir, "scevan_results2_complete.rds"))

# Print a confirmation message
cat("SCEVAN results saved to:", file.path(scevan_dir, "scevan_results2.txt"), "\n")

# save.image("analysis.scevan_19may.RData")

cat("✅ SCEVAN complete. Metadata added and session saved.\n")

# head(obj_scevan2@meta.data, 2)
#                                     orig.ident nCount_RNA nFeature_RNA
# LGG-04-1_LGG-04-1_AAACCCAAGTCCCAGC-1   LGG-04-1      10672         2982
# LGG-04-1_LGG-04-1_AAACCCAAGTGCCGAA-1   LGG-04-1       1781          874
#                                       sample percent.mt RNA_snn_res.0.5
# LGG-04-1_LGG-04-1_AAACCCAAGTCCCAGC-1 LGG-04-1  8.7706147               0
# LGG-04-1_LGG-04-1_AAACCCAAGTGCCGAA-1 LGG-04-1  0.9545199               0
#                                      seurat_clusters rpca_clusters  class
# LGG-04-1_LGG-04-1_AAACCCAAGTCCCAGC-1               0             0 normal
# LGG-04-1_LGG-04-1_AAACCCAAGTGCCGAA-1               8             8 normal
#                                     confidentNormal subclone
# LGG-04-1_LGG-04-1_AAACCCAAGTCCCAGC-1            <NA>       NA
# LGG-04-1_LGG-04-1_AAACCCAAGTGCCGAA-1            <NA>       NA
#
# colnames(scevan_results)
# [1] "class"           "confidentNormal" "subclone"  
# ----------------------------- #
#     Visualization & Figures   #
# ----------------------------- #
# table(obj_scevan2@meta.data$class)
# 
# filtered   normal    tumor 
#     467     6466     4122 
# table(obj_scevan2@meta.data$confidentNormal)

# yes 
# 30 
# table(obj_scevan2@meta.data$subclone)

#  1   2   3   4   5   6   7 
# 841 499 413 852 231 990 296

# orig.ident         # Sample ID
# nCount_RNA         # UMI count
# nFeature_RNA       # Number of genes detected
# sample             # Redundant with orig.ident
# percent.mt         # % mitochondrial reads
# RNA_snn_res.0.5    # Resolution-based clustering
# seurat_clusters    # Seurat default clusters
# rpca_clusters      # Clusters from RPCA integration
# class              # Likely CopyKAT class (e.g., normal, aneuploid)
# confidentNormal    # Possibly from CopyKAT/SCEVAN — not yet assigned
# subclone           # Possibly from SCEVAN — NA

table(obj_scevan2@meta.data$class)

# filtered   normal    tumor 
#      48      562      390 

######################################################################################
###################################################################################### display the data in
###################################################################################### obj_scevan2

# meta.data <- obj_scevan2@meta.data

setwd("../")
getwd()

# Ensure UMAP is computed
if (!"umap" %in% Reductions(obj_scevan2)) {
  obj_scevan2 <- RunUMAP(obj_scevan2, dims = 1:20)
}

# DimPlot() → for categorical data (e.g., cluster ID, cell type, condition)
# FeaturePlot() → for continuous data (e.g., CNA.score, percent.mt, nCount_RNA, gene expression)

# 1. UMAP by Seurat Clusters
p1 <- DimPlot(obj_scevan2, reduction = "umap", group.by = "seurat_clusters", label = TRUE,  pt.size = 0.5) +
      ggtitle("Seurat Clusters")  +
      theme(plot.title = element_text(size = 10)) 

ggsave(file = file.path(scevan_dir,"UMAP_seurat_clusters.png"), plot = p1, width = 6, height = 5, dpi = 300)

# 2. UMAP by Class
p2 <- DimPlot(obj_scevan2, reduction = "umap", group.by = "class", label = FALSE,  pt.size = 0.5) +
      ggtitle("Class")  +
      theme(plot.title = element_text(size = 10)) 

ggsave(file = file.path(scevan_dir,"UMAP_scevan_class.png"), plot = p2, width = 6, height = 5, dpi = 300)

# 3. UMAP by Confident Normal
p3 <- DimPlot(obj_scevan2, reduction = "umap", group.by = "confidentNormal", label = FALSE,  pt.size = 0.5) +
      ggtitle("Confident Normal")  +
      theme(plot.title = element_text(size = 10)) 

ggsave(file = file.path(scevan_dir,"UMAP_scevan_confidentNormal.png"), plot = p3, width = 6, height = 5, dpi = 300)

# 4. UMAP by SCEVAN Subclone
p4 <- DimPlot(obj_scevan2, reduction = "umap", group.by = "subclone", label = FALSE,  pt.size = 0.5) +
      ggtitle("SCEVAN Subclone")  +
      theme(plot.title = element_text(size = 10)) 

ggsave(file = file.path(scevan_dir,"UMAP_scevan_subclone.png"), plot = p4, width = 6, height = 5, dpi = 300)


# Combine p1 + p2
panel1 <- p1 | p2
ggsave(file = file.path(scevan_dir, "UMAP_scevan_cluster_class.png"), 
                        plot = panel1, width = 12, height = 5, dpi = 300)

# Combine p3 + p4
panel2 <- p3 | p4
ggsave(file = file.path(scevan_dir, "UMAP_scevan_confidentNormal_subclone.png"), 
                        plot = panel2, width = 12, height = 5, dpi = 300)

######################################################################################
######################################################################################
######################################################################################
###################################################################################### featurePlot
###################################################################################### CLONALITY 
######################################################################################
######################################################################################
######################################################################################

# Check if subclone information exists
if ("subclone" %in% colnames(obj_scevan2@meta.data)) {

  # Convert subclone factor to numeric
  obj_scevan2$subclone_numeric <- as.numeric(as.factor(obj_scevan2$subclone))

  # Create FeaturePlot
  p_subclone <- FeaturePlot(obj_scevan2,
                            features = "subclone",
                            cols = c("lightgrey", "purple"), pt.size = 0.5) +
    ggtitle("FeaturePlot: Subclone") + 
    theme(plot.title = element_text(size = 10)) 

  # Save the plot
  ggsave(file = file.path(scevan_dir,"scevan_subclone_featureplot.png"), 
                          plot = p_subclone, width = 6, height = 5, dpi = 300)

} else {
  warning("subclone column not found in meta.data")
}

######################################################################################
######################################################################################
###################################################################################### MATRICES
######################################################################################
# _CNAmtx.RData
# _CNAmtxSubclones.RData
# _count_mtx_annot.RData

# CNA_mtx_relat is the matrix of relative CNA values per cell, where:

# ✅ Rows → Genomic segments (bins or regions, e.g., chromosome arms, 100kb bins)#
# ✅ Columns → Single cells (or barcodes)
# ✅ Values → Relative copy number alterations (adjusted for normalization/reference)

setwd("./output")

######################################################################################
######################################################################################
###################################################################################### load the results

# above :

# load("_CNAmtx.RData")
# str(CNA_mtx_relat) load("_CNAmtxSubclones.RData")

# CNA_mtx_subclones = results.com 
# str(CNA_mtx_subclones) load("_count_mtx_annot.RData")

# str(count_mtx_annot) load("PlotOncoHeat.RData")
# str(oncoHeat)

CNA_mtx_clone  = CNA_mtx_subclones 

######################################################################################
###################################################################################### load the results
######################################################################################
###################################################################################### from virtual env

# enva <- new.env()
# load("_CNAmtx.RData", envir = enva)
# print(ls(enva))
# CNA_mtx_relat <- enva$CNA_mtx_relat

# head(colnames(CNA_mtx_relat),4)
# [1] "LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1"
# [2] "LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1"
# [3] "LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1"
# [4] "LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1"
# dim(CNA_mtx_relat)
# [1] 7166  952

# envb <- new.env()
# load("_CNAmtxSubclones.RData", envir = envb)
# print(ls(envb))

# CNA_mtx_clone = envb$results.com

# head(colnames(CNA_mtx_clone),4)
# [1] "LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1"
# [2] "LGG-04-1_LGG-04-1_AATAGAGAGACCTCAT-1"
# [3] "LGG-04-1_LGG-04-1_ACGGGTCCACCACATA-1"
# [4] "LGG-04-1_LGG-04-1_AGCTCAATCACCTCTG-1"
# dim(CNA_mtx_clone)
# [1] 7166  390

# envc <- new.env()
# load("_count_mtx_annot.RData", envir = envc)
# print(ls(envc))

# count_mtx_annot = envc$count_mtx_annot

# head(count_mtx_annot)
#                seqnames   start     end         gene_id gene_name
# ENSG00000188976        1  944204  959309 ENSG00000188976     NOC2L
# ENSG00000188290        1  998962 1000172 ENSG00000188290      HES4
# ENSG00000187608        1 1001138 1014541 ENSG00000187608     ISG15
# dim(count_mtx_annot)
# [1] 7166    5

# dim((CNA_mtx_relat))
# [1] 7166  952
# dim((results.com))
# [1] 7166  390
# dim((count_mtx_annot))
# [1] 7166    5

setwd("../")

######################################################################################
######################################################################################
######################################################################################
######################################################################################
###################################################################################### scevan_output_CNA_mtx_segm_mean_heatmap

# pdf(file = file.path(scevan_dir,"scevan_output_CNA_mtx_segm_mean_heatmap.pdf"), 
#    width = 10, height = 8)
# Heatmap(enva$CNA_mtx_relat, name = "segm_mean", show_row_names = FALSE, show_column_names = FALSE)
# Heatmap(CNA_mtx_relat, name = "segm_mean", show_row_names = FALSE, show_column_names = FALSE, cluster_rows = FALSE)
# Heatmap(CNA_mtx_relat, 
#        name = "segm_mean", 
#        show_row_names = FALSE, 
#        show_column_names = FALSE,
#        row_title = "Genome Segments",
#        column_title = "Cells",
#        row_title_side = "left",
#        column_title_side = "top",
#        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
#        column_title_gp = gpar(fontsize = 12, fontface = "bold"))
# dev.off()

library(ComplexHeatmap)

# Build the heatmap object once
ht <- Heatmap(
  CNA_mtx_relat, 
  name = "segm_mean", 
  show_row_names = FALSE, 
  show_column_names = FALSE,
  row_title = "Genome Segments",
  column_title = "Cells",
  row_title_side = "left",
  column_title_side = "top",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

# Save to PDF
pdf(file = file.path(scevan_dir, "scevan_output_CNA_mtx_segm_mean_heatmap.pdf"), width = 10, height = 8)
draw(ht)
dev.off()

# Save to PNG
png(file = file.path(scevan_dir, "scevan_output_CNA_mtx_segm_mean_heatmap.png"), width = 1600, height = 1200, res = 200)
draw(ht)
dev.off()

######################################################################################
######################################################################################
######################################################################################
######################################################################################
###################################################################################### scevan_output_CNA_mtx_segm_mean_heatmap_subclones

# pdf(file = file.path(scevan_dir,"scevan_output_CNA_mtx_segm_mean_heatmap_subclones.pdf"), 
#    width = 10, height = 8)
# Heatmap(envb$results.com, name = "segm.mean", show_row_names = FALSE, show_column_names = FALSE)
# Heatmap(CNA_mtx_clone, name = "segm.mean", show_row_names = FALSE, show_column_names = FALSE, cluster_rows = FALSE)
# Heatmap(CNA_mtx_clone, 
#        name = "segm.mean", 
#        show_row_names = FALSE, 
#        show_column_names = FALSE,
#        row_title = "Genome Segments",
#        column_title = "Cells",
#        row_title_side = "left",
#        column_title_side = "top",
#        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
#        column_title_gp = gpar(fontsize = 12, fontface = "bold"))
# dev.off()

library(ComplexHeatmap)

# Build the heatmap object once
ht_subclone <- Heatmap(
  CNA_mtx_clone, 
  name = "segm.mean", 
  show_row_names = FALSE, 
  show_column_names = FALSE,
  row_title = "Genome Segments",
  column_title = "Cells",
  row_title_side = "left",
  column_title_side = "top",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

# Save to PDF
pdf(file = file.path(scevan_dir, "scevan_output_CNA_mtx_segm_mean_heatmap_subclones.pdf"), width = 10, height = 8)
draw(ht_subclone)
dev.off()

# Save to PNG
png(file = file.path(scevan_dir, "scevan_output_CNA_mtx_segm_mean_heatmap_subclones.png"), width = 1600, height = 1200, res = 200)
draw(ht_subclone)
dev.off()

######################################################################################
######################################################################################
######################################################################################
######################################################################################

# You’re plotting a histogram of all segment means across all cells.
# CNA_mtx_relat in SCEVAN (or enva$CNA_mtx_relat) is a matrix of segment means for each cell across many genomic segments. 
# So:
# Rows = cells
# Columns = segments
# If you have 1,000 cells and 1,000 segments, that matrix has 1,000,000 (1 million) values.
# The frequencies in the histogram are counts of these 1 million values, not the number of cells.

# Plot histogram
# pdf(file = file.path(scevan_dir,"scevan_output_CNA_mtx_segm_mean_histogram.pdf"), 
#    width = 7, height = 5)
# ggplot(data.frame(segm_mean = as.vector(enva$CNA_mtx_relat)), 
# ggplot(data.frame(segm_mean = as.vector(CNA_mtx_relat)), 
#  aes(x = segm_mean)) +
#  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
#  theme_minimal() +
#  labs(title = "CNA segm.mean values per genome bin",
#       x = "segm_mean",
#       y = "Frequency")
# dev.off()

# Plot histogram
hist_plot <- ggplot(data.frame(segm_mean = as.vector(CNA_mtx_relat)), 
                    aes(x = segm_mean)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "CNA segm.mean values per genome bin",
    x = "segm_mean",
    y = "Frequency"
  )

# Save to PDF
pdf(file = file.path(scevan_dir, "scevan_output_CNA_mtx_segm_mean_histogram.pdf"), width = 7, height = 5)
print(hist_plot)
dev.off()

# Save to PNG
png(file = file.path(scevan_dir, "scevan_output_CNA_mtx_segm_mean_histogram.png"), width = 1200, height = 900, res = 200)
print(hist_plot)
dev.off()

# Plot histogram
# pdf(file = file.path(scevan_dir,"scevan_output_CNA_mtx_segm_mean_histogram_subclones.pdf"), 
#    width = 7, height = 5)

# ggplot(as.vector(envb$results.com), 
# ggplot(data.frame(segm_mean = as.vector(CNA_mtx_clone)), 
#  aes(x = segm_mean)) +
#  geom_histogram(bins = 100, fill = "skyblue", color = "black") +
#  theme_minimal() +
#  labs(title = "Subclonal CNA segm.mean values per genome bin",
#       x = "segm.mean",
#       y = "Frequency")
# dev.off()

library(ggplot2)

# Create the ggplot object
hist_subclone_plot <- ggplot(data.frame(segm_mean = as.vector(CNA_mtx_clone)), 
                             aes(x = segm_mean)) +
  geom_histogram(bins = 100, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Subclonal CNA segm.mean values per genome bin",
    x = "segm.mean",
    y = "Frequency"
  )

# Save to PDF
pdf(file = file.path(scevan_dir, "scevan_output_CNA_mtx_segm_mean_histogram_subclones.pdf"), width = 7, height = 5)
print(hist_subclone_plot)
dev.off()

# Save to PNG
png(file = file.path(scevan_dir, "scevan_output_CNA_mtx_segm_mean_histogram_subclones.png"), width = 1200, height = 900, res = 200)
print(hist_subclone_plot)
dev.off()

###################################################################################### per each CLONE
######################################################################################
###################################################################################### the histogram of seg mean and CN
###################################################################################### for each clone
######################################################################################

setwd("./output")

# files <- c(
#  "_Clonal_CN.seg",
#  "_subclone1_CN.seg",
#  "_subclone2_CN.seg",
#  "_subclone3_CN.seg",
#  "_subclone4_CN.seg",
#  "_subclone5_CN.seg"
# )

# clonal segmentation file : _Clonal_CN.seg

#  Chr       Pos       End CN segm.mean
# 1   1    944204 154220628  2 -0.127514
# 2   1 154220179 248849517  2  0.015655
# 3   2    217730 241817413  2  0.024638
# 4   3    196596  58295693  2  0.016539
# 5   3  58306247 119241103  4  0.141726
# 6   3 119294373 170860380  2 -0.03005

# Create the combined segmentation dataframe from the input list
seg_list <- lapply(names(seg_files_data), function(clone_name) {
  df <- seg_files_data[[clone_name]]
  colnames(df) <- c("chr", "start", "end", "CN", "segm.mean")
  df$clone <- clone_name
  df
})

# Bind all clones together
seg_all <- bind_rows(seg_list)

# Loop over clones and generate plots
unique_clones <- unique(seg_all$clone)

for (clone in unique_clones) {
  df <- seg_all %>% filter(clone == !!clone)

  # Histogram of segm.mean
  p1 <- ggplot(df, aes(x = segm.mean)) +
    geom_histogram(bins = 100, fill = "#9B59B6", color = "black") +
    scale_x_continuous(
      breaks = seq(floor(min(df$segm.mean, na.rm = TRUE)), ceiling(max(df$segm.mean, na.rm = TRUE)), by = 0.1),
      expand = c(0, 0.05)
    ) +
    theme_minimal(base_size = 12) +
    labs(title = "segm.mean", x = "segm.mean", y = "Frequency")

  # Histogram of CN
  p2 <- ggplot(df, aes(x = CN)) +
    geom_histogram(binwidth = 1, fill = "#FAA43A", color = "black", center = 0.1) +
    scale_x_continuous(
      breaks = seq(floor(min(df$CN, na.rm = TRUE)), ceiling(max(df$CN, na.rm = TRUE)), by = 0.5),
      expand = c(0, 0.05)
    ) +
    theme_minimal(base_size = 12) +
    labs(title = "Copy Number", x = "Copy Number", y = "Frequency")

  # Combine plots side by side
  final_plot <- p1 + p2 +
    plot_layout(ncol = 2, widths = c(0.8, 0.8)) +
    plot_annotation(
      title = paste("Segmentation Histograms: Clone", clone),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )

  # Save to PDF
  ggsave(filename = paste0(clone, "_segmentation_histograms.pdf"),
         plot = final_plot, width = 10, height = 5)

  # Save to PNG
  ggsave(filename = paste0(clone, "_segmentation_histograms.png"),
         plot = final_plot, width = 10, height = 5, dpi = 300)
}

######################################################################################
######################################################################################
###################################################################################### How can we use ? 
######################################################################################
######################################################################################

# segmentation file, 
# matrix file, and 
# metadata file

# The CNA matrix (CNA_mtx_relat) has:

# Rows = genomic bins or segments
# Columns = single cells
# Values = relative CNV signals (like log2 ratios)

dim((CNA_mtx_relat))
# [1] 7166  952
dim((CNA_mtx_clone))
# [1] 7166  390
dim((count_mtx_annot))
# [1] 7166    5

# To note that these data structures have the same number of rows (7166 in this example).

# 🧬 CNA_mtx_relat
# In SCEVAN, CNA_mtx_relat typically contains the relative copy number alteration (CNA) values 
# for each segment (rows) and each single cell (columns).
# It’s a numeric matrix with no explicit gene annotations—just values of CNAs.

# 🗂️ count_mtx_annot
# count_mtx_annot is a matrix that usually includes annotations for each bin/segment:
# Chromosome
# Start and End coordinates
# Gene names (if available)

# This matrix is often used as input for computing CNA_mtx_relat (info by GPT4 :) 

# ✅ So:
# Annotations are in count_mtx_annot
# CNA values are in CNA_mtx_relat
# Both can be matched by their row order (if no changes were made to row ordering).

# obj_scevan2
# An object of class Seurat 
# 36601 features across 1000 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 7 layers present: data, counts, scale.data.LGG-04-1, scale.data.LGG-04-2, scale.data.LGG-04-3, scale.data.LGG-03, scale.data
# 4 dimensional reductions calculated: pca, harmony, umap.harmony, umap

CNA_mtx_relat_annot = cbind(count_mtx_annot, CNA_mtx_relat)
CNA_mtx_clone_annot = cbind(count_mtx_annot, CNA_mtx_clone)

CNA_mtx_relat_annot[1:4, 1:5]
CNA_mtx_clone_annot[1:4, 1:5]

# CNA_mtx_relat_annot[1:4, 1:7]
#                seqnames   start     end         gene_id gene_name
# ENSG00000188976        1  944204  959309 ENSG00000188976     NOC2L
# ENSG00000188290        1  998962 1000172 ENSG00000188290      HES4
# ENSG00000187608        1 1001138 1014541 ENSG00000187608     ISG15
# ENSG00000078808        1 1216908 1232031 ENSG00000078808      SDF4
#                 LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1
# ENSG00000188976                           -0.1706296
# ENSG00000188290                           -0.1706296
# ENSG00000187608                           -0.1706296
# ENSG00000078808                           -0.1706296
#                 LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1
# ENSG00000188976                          -0.02180575
# ENSG00000188290                          -0.02180575
# ENSG00000187608                          -0.02180575
# ENSG00000078808                          -0.02180575


# CNA_mtx_clone_annot[1:4, 1:7]
#                seqnames   start     end         gene_id gene_name
# ENSG00000188976        1  944204  959309 ENSG00000188976     NOC2L
# ENSG00000188290        1  998962 1000172 ENSG00000188290      HES4
# ENSG00000187608        1 1001138 1014541 ENSG00000187608     ISG15
# ENSG00000078808        1 1216908 1232031 ENSG00000078808      SDF4
#                LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1
# ENSG00000188976                           -0.1698048
# ENSG00000188290                           -0.1698048
# ENSG00000187608                           -0.1698048
# ENSG00000078808                           -0.1698048
#                LGG-04-1_LGG-04-1_AATGGAAGTAACTAAG-1
# ENSG00000188976                           -0.3292702
# ENSG00000188290                           -0.3292702
# ENSG00000187608                           -0.3292702
# ENSG00000078808                           -0.3292702

 head(obj_scevan2@meta.data, 5)
#                                     orig.ident nCount_RNA nFeature_RNA
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1   LGG-04-1        517          265
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1   LGG-04-1       7999         2722
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1   LGG-04-1       9911         3005
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1   LGG-04-1      16196         3761
# LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1   LGG-04-1       8531         3480
#                                       sample percent.mt RNA_snn_res.0.5
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1 LGG-04-1   6.576402               8
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1 LGG-04-1   8.738592               0
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1 LGG-04-1   8.858844               0
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1 LGG-04-1   6.693011               2
# LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1 LGG-04-1   9.858164               1
#                                     seurat_clusters harmony_clusters    class
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1              14               14 filtered
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1               0                0   normal
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1               0                0   normal
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1               0                0   normal
# LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1               1                1    tumor
#                                     confidentNormal subclone subclone_numeric
# LGG-04-1_LGG-04-1_AAAGGTATCACTGTCC-1            <NA>       NA               NA
# LGG-04-1_LGG-04-1_AACGTCACAACGGCCT-1            <NA>       NA               NA
# LGG-04-1_LGG-04-1_AAGCATCTCCTCGATC-1            <NA>       NA               NA
# LGG-04-1_LGG-04-1_AAGGAATCACGCCACA-1            <NA>       NA               NA
# LGG-04-1_LGG-04-1_AAGTACCTCTCAACCC-1            <NA>        1                1

metadata = obj_scevan2@meta.data

######################################################################################
######################################################################################
######################################################################################
###################################################################################### CNV profile for one cell 

#  Line plot of CNV profile (for one cell)

cell_idx <- 101  # choose cell index
cell_cnv <- CNA_mtx_relat[, cell_idx]

# Save plot to PNG file
png(filename = paste0("cell", cell_idx, "_CNV_profile.png"),
    width = 1200, height = 400)

# Plot with thicker line
plot(cell_cnv, type = "l",
     col = "blue",
     lwd = 2,           # line width for thicker line
     xlab = "Genome Bin",
     ylab = "CNV Signal",
     main = paste("CNV Profile for Cell", cell_idx))
abline(h = 0, lty = 2)

# Close the file
dev.off()

######################################################################################
######################################################################################
######################################################################################
###################################################################################### CNA_mtx_segm_mean_heatmap_clusters

# heatmap of CNA data from single-cell RNA-seq, showing:
# CNAs across segments (rows) for each cell (columns).
# Cells annotated by sample, class, and cluster.
# Colored by segm.mean values (CNA values).

# 1. Load CNA matrix and metadata
cna_matrix <- CNA_mtx_relat  # matrix: rows = segments, columns = cells
meta <- obj_scevan2@meta.data  # Seurat metadata

# 2. Remove cells with class == "filtered"
meta <- meta[meta$class != "filtered", , drop = FALSE]

# 3. Subset and align the matrix with filtered metadata
cna_matrix <- cna_matrix[, colnames(cna_matrix) %in% rownames(meta)]
meta <- meta[colnames(cna_matrix), , drop = FALSE]  # align meta order to matrix

# 4. Convert harmony_clusters to factor for categorical coloring
meta$harmony_clusters <- as.factor(meta$harmony_clusters)

# 5. Generate discrete colors for sample and cluster annotations
sample_levels <- unique(meta$sample)
sample_colors <- setNames(
  brewer.pal(n = max(3, length(sample_levels)), name = "Set3")[seq_along(sample_levels)],
  sample_levels
)

cluster_levels <- levels(meta$harmony_clusters)
cluster_colors <- setNames(
  brewer.pal(n = max(3, length(cluster_levels)), name = "Paired")[seq_along(cluster_levels)],
  cluster_levels
)

# 6. Create top column annotation
column_annot <- HeatmapAnnotation(
  sample = meta$sample,
  class = meta$class,
  cluster = meta$harmony_clusters,
  col = list(
    sample = sample_colors,
    class = c("normal" = "#1b9e77", "tumor" = "#7570b3"),
    cluster = cluster_colors
  ),
  show_annotation_name = TRUE
)

# 7. Define custom color gradient for segm.mean values
col_fun <- colorRamp2(
  breaks = c(-1, -0.4, -0.2, 0, 0.2, 0.4, 1),
  colors = c("#2166AC", "#67A9CF", "#C6DBEF", "#F7F7F7", "#FDAE6B", "#EF8A62", "#B2182B")
)

# 8. Draw the heatmap
# Create the heatmap object
ht <- Heatmap(
  cna_matrix,
  name = "segm.mean",
  col = col_fun,
  top_annotation = column_annot,
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  column_title = "Cell annotations",
  row_title = "Genome Segments",
  column_title_side = "top",
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

# Save to PDF
pdf("scevan_output_CNA_mtx_segm_mean_heatmap_clusters.pdf", width = 12, height = 8)
draw(ht)
dev.off()

# Save to PNG
png("scevan_output_CNA_mtx_segm_mean_heatmap_clusters.png", width = 1600, height = 1200, res = 200)
draw(ht)
dev.off()

######################################################################################
######################################################################################
######################################################################################
###################################################################################### CNA_mtx_segm_mean_heatmap_clones

# Clustered heatmap that highlights:
# ✅ CNA values (segm.mean) for each segment (rows) across cells (columns)
# ✅ Sample, class (normal/tumor), and cluster metadata annotations for each cell
# ✅ Clonal structure (since you’re using CNA_mtx_clone, likely indicating clone-based CNA matrix)
# ✅ The data is hierarchically clustered for pattern discovery

# 1. Load CNA matrix and metadata
cna_matrix <- CNA_mtx_clone    # matrix: rows = segments, columns = cells
meta <- obj_scevan2@meta.data  # Seurat metadata

# 2. Remove cells with class == "filtered"
meta <- meta[meta$class != "filtered", , drop = FALSE]

# 3. Subset and align the matrix with filtered metadata
cna_matrix <- cna_matrix[, colnames(cna_matrix) %in% rownames(meta)]
meta <- meta[colnames(cna_matrix), , drop = FALSE]  # align meta order to matrix

# 4. Convert harmony_clusters to factor for categorical coloring
meta$harmony_clusters <- as.factor(meta$harmony_clusters)

# 5. Generate discrete colors for sample and cluster annotations
sample_levels <- unique(meta$sample)
sample_colors <- setNames(
  brewer.pal(n = max(3, length(sample_levels)), name = "Set3")[seq_along(sample_levels)],
  sample_levels
)

cluster_levels <- levels(meta$harmony_clusters)
cluster_colors <- setNames(
  brewer.pal(n = max(3, length(cluster_levels)), name = "Paired")[seq_along(cluster_levels)],
  cluster_levels
)

# 6. Create top column annotation
column_annot <- HeatmapAnnotation(
  sample = meta$sample,
  class = meta$class,
  cluster = meta$harmony_clusters,
  col = list(
    sample = sample_colors,
    class = c("normal" = "#1b9e77", "tumor" = "#7570b3"),
    cluster = cluster_colors
  ),
  show_annotation_name = TRUE
)

# 7. Define custom color gradient for segm.mean values
col_fun <- colorRamp2(
  breaks = c(-1, -0.4, -0.2, 0, 0.2, 0.4, 1),
  colors = c("#2166AC", "#67A9CF", "#C6DBEF", "#F7F7F7", "#FDAE6B", "#EF8A62", "#B2182B")
)

# 8. Draw the heatmap
ht_clones <- Heatmap(
  cna_matrix,
  name = "segm.mean",
  col = col_fun,
  top_annotation = column_annot,
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  column_title = "Cell Annotations",
  row_title = "Genome Segments",
  column_title_side = "top",
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

# Save to PDF
pdf("scevan_output_CNA_mtx_segm_mean_heatmap_clones.pdf", width = 12, height = 8)
draw(ht_clones)
dev.off()

# Save to PNG
png("scevan_output_CNA_mtx_segm_mean_heatmap_clones.png", width = 1600, height = 1200, res = 200)
draw(ht_clones)
dev.off()

######################################################################################
######################################################################################
###################################################################################### chromosome-wise CNA heatmaps for each distance metric.
######################################################################################
######################################################################################
######################################################################################

# Extracts CNA data for that chromosome’s segments.

#  Creates a heatmap (Heatmap):

# ✅ CNA matrix for that chromosome
# ✅ Clustering for rows (cells), no clustering for columns (segments in order)
# ✅ Hides row and column names for clarity
# ✅ Creates chromosome-wise CNA heatmaps for each distance metric.
# ✅ Clusters cells based on CNA profiles using different distance metrics.
# ✅ Adds cell-level annotations (class and cluster).

# 📊 Rows: Cells
# Each row in the heatmap represents a single cell from the dataset.
# The ordering of rows is determined by clustering of the CNA profiles (e.g., hierarchical clustering based on Euclidean or cosine distances).

# 📊 Columns: Genomic Segments (per chromosome)
# Each column represents a genomic segment (like a bin/region of the genome, e.g., from CNA segmentation).

# The segments are grouped and ordered by chromosome:
# All segments from chromosome 1 first, then chromosome 2, and so on (like a genome-wide tiling).
# No clustering on columns — they’re ordered as they appear in the segmentation data.

# 🔎 Details for clarity
# Per-chromosome heatmaps:
# Each chromosome’s segments are visualized separately (separate smaller heatmaps), but combined into a single composite figure.

# Color gradient (segm.mean values):
# The heatmap’s colors show the segmentation mean values (CNA values) for each cell across each segment.

# Row annotations:
# Color bars on the left of the heatmap indicate cell-level class (like normal/tumor) and harmony_clusters membership.

# Input data: assumed already in workspace

# CNA_mtx_relat: matrix of segm.mean (rows = segments/genes, columns = cells)
# count_mtx_annot: data frame with chromosome, start, end, gene_id, gene_name

library(RColorBrewer)
library(circlize)
library(dplyr)
library(fastcluster)
library(parallelDist)
library(pheatmap)

# Capability of pheatmap:

# Usually fine up to about 10,000 rows and 2,000 columns.
# Beyond that, it can become very slow or even fail (depending on your system RAM).

# Capability of ComplexHeatmap:
# it’s in the upper range of typical usage

# The key libraries that can significantly speed up large heatmap creation are:
# 1. fastcluster - Much faster hierarchical clustering (up to 10x faster)
# 2. parallelDist - Parallel distance computation using multiple CPU cores
# 3. pheatmap - Often much faster than ComplexHeatmap for large matrices

# Speed improvements:

# fastcluster: Optimized C++ implementation of hierarchical clustering
# parallelDist: Uses all available CPU cores for distance calculations
# pheatmap: Streamlined implementation that's often 2-5x faster than ComplexHeatmap for large datasets
# Cairo graphics: Often faster PNG rendering

# Available Distance Metrics in parallelDist::parDist()
# --------------------------
# "euclidean" – default
# "manhattan" – city block / taxicab
# "maximum" – Chebyshev / L∞
# "canberra" – more sensitive to small differences
# "binary" – for binary data
# "minkowski" – generalized distance (requires p parameter)
# "cosine" – cosine distance
# "correlation" – correlation-based distance
# "hamming" – Hamming distance (binary)
# "jaccard" – Jaccard dissimilarity (binary)
# "spearman" – Spearman rank correlation distance
# "kendall" – Kendall’s tau correlation distance

# Available Linkage Methods:
# --------------------------
# "ward.D"  Ward’s minimum variance method (minimizes total within-cluster variance).
# "ward.D2" Like "ward.D" but uses squared Euclidean distances.
# "single"  Nearest point linkage (single linkage, minimum distance).
# "complete"  Farthest point linkage (maximum distance).
# "average" Average linkage (UPGMA — unweighted pair-group method with arithmetic mean).
# "mcquitty"  McQuitty’s method (WPGMA — weighted pair-group method with arithmetic mean).
# "median"  Median linkage (WPGMC — weighted pair-group method using centroids).
# "centroid"  Centroid linkage (UPGMC — unweighted pair-group method using centroids).

# metadata - your cell metadata
metadata = obj_scevan2@meta.data

# Extract genomic annotation columns (first 5 columns)
genomic_info <- CNA_mtx_relat_annot[, 1:5]

# Extract CNA values (columns 6 onwards - the cell data)
cna_matrix <- CNA_mtx_relat_annot[, 6:ncol(CNA_mtx_relat_annot)]

######################################################################################
######################################################################################
###################################################################################### display the HEATMAP of Segm Mean
######################################################################################
###################################################################################### code generated by Claude Sonnet
######################################################################################
###################################################################################### tumor and normal

source("../script.display.genome_wide_segm_means_heatmaps_per_chromosome.normal.and.tumor.R")
source("../script.display.genome_wide_segm_means_heatmaps_per_chromosome.only.tumor.R")

######################################################################################
######################################################################################
###################################################################################### Effective methods 
###################################################################################### to reduce
###################################################################################### the Matrix Size
###################################################################################### 

# There are 5 effective methods to reduce your matrix size while preserving biological information:

# 1. Genomic Binning (Recommended for CNA)

# What: Groups adjacent genes into larger genomic windows (e.g., 1MB bins)
# Preserves: Spatial chromosome structure and large-scale CNAs
# Reduction: ~7166 genes → ~3000 bins (depending on bin size)
# Best for: Detecting chromosomal-level alterations

# 2. Principal Component Analysis (PCA)

# What: Transforms genes into principal components that capture variance
# Preserves: 95% of variance with ~100-200 components
# Reduction: ~7166 genes → ~100-200 PCs
# Best for: Preserving cell-to-cell differences

# 3. Variance-based Filtering

# What: Keeps only genes with highest variance across cells
# Preserves: Most informative genes for distinguishing cell types
# Reduction: ~7166 genes → ~2000 most variable genes
# Best for: Focusing on genes that differ between cells

# 4. Segmentation-based Reduction
#
# What: Merges adjacent genes with similar CNA patterns
# Preserves: Natural biological segments/breakpoints
# Reduction: Variable, depends on your data's segmentation
# Best for: Maintaining biological CNV boundaries

# 5. Cytoband-based Reduction

# What: Groups genes by chromosomal arms/bands (p and q arms)
# Preserves: Cytogenetic structure
# Reduction: ~7166 genes → ~200-300 cytobands
# Best for: Clinical/cytogenetic interpretation

# Bin Size :

# 0.1 MB Gene-level CNAs
# 0.5 MB Focal CNAs
# 1 MB Balanced
# 2 MB Arm-level CNAs
# 5 MB Quick preview

# Summary: Computing Segment Mean for Cytobands
# When binning genes into cytobands, you aggregate all gene-level segment means within each cytoband region. Here are the main approaches:

# 1. Simple Mean (Most Common) ⭐
# For each cytoband, average all gene segment means
# cytoband_segmean = mean(gene1_segmean, gene2_segmean, ..., geneN_segmean)
# Best for: Balanced gene distribution, standard analysis

# 2. Weighted Mean by Gene Length
# Weight by how much genomic space each gene covers
# cytoband_segmean = weighted.mean(gene_segmeans, weights = gene_lengths)
# Best for: Regions with very different gene sizes

# 3. Median (Robust)
# Use median instead of mean to reduce outlier influence
# cytoband_segmean = median(gene1_segmean, gene2_segmean, ..., geneN_segmean)
# Best for: Noisy data with potential outlier genes

# 4. Coverage-Weighted Mean
# Weight by how much of each gene actually overlaps the cytoband
# cytoband_segmean = weighted.mean(gene_segmeans, weights = overlap_lengths)
# Best for: Genes that partially overlap cytoband boundaries
# My Recommendation

# For your SCEVAN data, I recommend the simple mean approach because:

# Standard practice in CNA analysis
# Computationally efficient
# Interpretable results
# Works well when genes are reasonably distributed

# Total Number of Cytobands in Human Genome
# The human genome contains approximately 862 cytobands when subdivided into the finest resolution sub-bands UcscUcsc.

# Here's the breakdown:
# Resolution Levels:

# ~500 Major Bands: At standard resolution, there are about 500 bands visible across the entire genome UcscUcsc
# 862 Sub-bands: A final subdivision into a total of 862 sub-bands is made by adding a period and another digit to the band, 
# resulting in naming like 3p26.3, 3p26.2, etc. 

# Cytoband Naming System:

# Chromosome: 1-22, X, Y (24 chromosomes total)
# Arms: p (short arm) and q (long arm)
# Major bands: numbered sequentially from centromere outward
# Sub-bands: decimal subdivisions (e.g., 3p26.1, 3p26.2, 3p26.3)

# For Your CNA Analysis:
# If you use cytoband-based binning instead of 1MB bins:

# Reduction: 7166 genes → ~862 cytobands (8.3x reduction)
# Speed: Much faster clustering than 1MB bins
# Biology: Clinically relevant cytogenetic resolution
# Interpretation: Familiar to clinical geneticists

# Comparison with 1MB binning:

# 1MB bins: ~1600 bins (4.5x reduction)
# Cytobands: ~862 bands (8.3x reduction)

# 640 bins of 5 MB ; 320 bins of 10 MB

######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
###################################################################################### cytobands matrix

# Step 1 : Download UCSC Cytobands for hg38

# Simple function to download cytoband data
# cat("Downloading UCSC cytoband data for hg38...\n")

# UCSC URL for hg38 cytobands
# cytoband_url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"

# Download to a temporary file
# temp_file <- tempfile(fileext = ".txt.gz")
# download.file(cytoband_url, temp_file, mode = "wb", quiet = TRUE)

# Read the cytoband data
# cytoband_data <- read.table(temp_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Set proper column names (this is the UCSC format)
# colnames(cytoband_data) <- c("seqnames", "start", "end", "cytoband_name", "gieStain")

# Remove 'chr' prefix to match your data format (1, 2, 3 instead of chr1, chr2, chr3)
# cytoband_data$seqnames <- gsub("^chr", "", cytoband_data$seqnames)

# Keep only chromosomes 1-22 (remove X, Y, and any others)
# cytoband_data <- cytoband_data[cytoband_data$seqnames %in% as.character(1:22), ]

# Clean up the temporary file
# unlink(temp_file)

# Show what we downloaded
# cat("SUCCESS! Downloaded", nrow(cytoband_data), "cytobands for chromosomes 1-22\n")
# cat("First few cytobands:\n")
# print(head(cytoband_data))

# cat("\nChromosomes included:\n")
# print(sort(as.numeric(unique(cytoband_data$seqnames))))

# cat("\nCytoband data is now stored in: cytoband_data\n")
# cat("Ready for next step!\n")

# Step 2: Match Genes to Cytobands and Calculate Robust Averages
# Step 3: Computing a cytoband matrix to visualize on heatmaps  

# head(cytoband_data,3)
#  seqnames   start     end cytoband_name gieStain
# 1        1       0 2300000        p36.33     gneg
# 2        1 2300000 5300000        p36.32   gpos25
# 3        1 5300000 7100000        p36.31     gneg

# Common Reasons for Fewer Cytobands:

# 1. Empty Cytobands (0 genes)
# Your SCEVAN data doesn't cover all genome regions
# Some cytobands are in gene-poor areas (heterochromatin)
# SCEVAN focuses on gene-rich regions
# 2. Small Cytobands (<3 genes)
# Filtered out for statistical reliability
# Can't reliably average just 1-2 genes
# Current threshold: minimum 3 genes per cytoband
# 3. Gene Density Variation
# Gene distribution is uneven across genome
# Some cytobands are naturally gene-dense
# Others are gene-poor (normal!)

######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
################################################################# 
################################################################# to choose a specific region from a seg file :
################################################################# 

# head(count_mtx_annot)
#                seqnames   start     end         gene_id gene_name
# ENSG00000188976        1  944204  959309 ENSG00000188976     NOC2L
# ENSG00000188290        1  998962 1000172 ENSG00000188290      HES4
# ENSG00000187608        1 1001138 1014541 ENSG00000187608     ISG15

# https://github.com/AntonioDeFalco/SCEVAN
# If you want to plot CN information at the single-cell level, you can obtain the region of the alteration of interest 
# from the *.seg file and plot the inferred CN ratio from CNA matrix, for example, like this:
# is computing mean CNV values for a specific region on chromosome 3 across all single cells.

# Why do we get the average ?
#
# In scRNA-seq CNV inference, especially with tools like SCEVAN or CopyKAT, you typically don’t get one single seg.mean 
# per cell per region. Instead, the CNA matrix (CNA_mtx_relat) is often:
# Cells × segments (or bins or genes)
# With many small segments per chromosome arm (e.g., 10–100+ bins just for chr3q)

# So if you're interested in a specific region (e.g., chr3q29), you don’t have one value per cell for that region 
# — you have multiple segment values per cell covering that region.

# If you have a segmentation-level summary per cell (like from CopyKAT SEG files or bulk GISTIC), 
# then yes — you might already have a seg.mean for chr3q29 per cell or clone, and averaging wouldn’t be needed. 
# But in most scRNA-seq CNA matrices, averaging over bins is essential.

######################################################################################
######################################################################################
###################################################################################### display a genome segment on chr 3
###################################################################################### Segm Mean

# Why multiple values per region per cell?

# 1. Segmented inference across the genome
# scRNA-seq CNV tools (e.g., SCEVAN, InferCNV, CopyKAT) infer CNVs by analyzing expression patterns across the genome.

# Instead of predefining large cytobands (e.g., chr3q29), they:
# Divide the genome into many smaller regions (segments or bins) — 
# often based on gene windows or fixed genomic intervals (e.g., 100kb–1Mb).

# Infer copy number per segment for each individual cell.

# As a result, for chr3q29, which spans ~36 Mb (from ~158 Mb to ~194 Mb), you’ll typically have:
# ~10–50+ bins or segments per cell, each with its own log-ratio or CNA score

# 2. scRNA-seq is noisy and sparse
# Due to dropout and low coverage, you can’t confidently call CNVs at the gene or small region level in each cell.
# Averaging across multiple segments reduces noise and yields a more robust summary per region per cell.

# Calculate mean CNA for chromosome 3 region
chr3reg <- apply(CNA_mtx_relat[count_mtx_annot$seqnames == 3 & 
                              count_mtx_annot$start >= 158644278 & 
                              count_mtx_annot$end <= 194498364,], 2, mean)

# Optional: View the selected genes (uncomment if you want to see them)
# filtered_genes <- count_mtx_annot$seqnames == 3 & 
#                   count_mtx_annot$start >= 158644278 & 
#                   count_mtx_annot$end <= 194498364
# selected_genes <- count_mtx_annot[filtered_genes, ]
# cat("Number of genes selected:", sum(filtered_genes), "\n")
# print(head(selected_genes))
# # Applies the mean function column-wise (because 2 = margin for columns)
# Match cell order and prepare for Seurat
chr3reg <- chr3reg[rownames(obj_scevan2@meta.data)]
names(chr3reg) <- rownames(obj_scevan2@meta.data)

# Convert to data frame with proper column name
chr3reg_df <- data.frame(chr3reg = chr3reg)
rownames(chr3reg_df) <- names(chr3reg)

# Add to Seurat metadata
obj_scevan2 <- AddMetaData(obj_scevan2, metadata = chr3reg_df)

# Create feature plot
# Seurat::FeaturePlot(obj_scevan2, "chr3reg", cols = c("gray", "red"))

# Save plot
png("example_UMAP_chr3_region_CNA_FeaturePlot_Segmn_Mean.png", width = 800, height = 800, res = 150)
Seurat::FeaturePlot(obj_scevan2, "chr3reg", cols = c("gray", "red"), pt.size = 0.3) +
  theme(
    plot.title = element_text(size = 8),        # Title font size
    axis.title = element_text(size = 6),        # Axis titles (UMAP_1, UMAP_2)
    axis.text = element_text(size = 6),         # Axis tick labels (numbers)
    legend.title = element_text(size = 6),      # Legend title font size
    legend.text = element_text(size = 5)        # Legend labels font size
  ) +
  labs(color = "Segmn Mean")                    # Custom legend title
dev.off()

# Excellent — your code is showing a UMAP heatmap of average CNA signal (segm.mean) 
# over a specific region on chromosome 3, for each single cell in your Seurat object.

# 🔍 What You're Showing:
# You're displaying, per cell, the:

# Mean CNA value (segm.mean) across all genomic segments in
# chr3:158,644,278–194,498,364
# (a ~36 Mb region on chromosome 3)

# 🔬 Biological Interpretation:
# Positive values (red) = copy number gain in that region
# Negative values (gray) = copy number loss or diploid

# Color gradient (gray → red) reflects the degree of amplification in that chr3 region
# Cells with higher mean CNA values in that region (i.e., amplified chr3) are shown in red on the UMAP.
# This plots the chr3 CNV value on a UMAP or tSNE, using a color gradient.
# Low signal (deletions) will appear gray, and high signal (amplifications) will appear red.

######################################################################################
######################################################################################
###################################################################################### display a genome segment on chr 3
###################################################################################### Gain / Losses

# Classify as gain/loss/neutral based on thresholds
obj_scevan2$chr3_cna <- ifelse(obj_scevan2$chr3reg > 0.2, "Gain",
                        ifelse(obj_scevan2$chr3reg < -0.2, "Loss", "Neutral"))

# Visualize on UMAP
png("example_UMAP_chr3_region_CNA_DimPlot_gain_loss_neutral.png", width = 800, height = 800, res = 150)
Seurat::DimPlot(obj_scevan2, group.by = "chr3_cna", cols = c("blue", "gray", "red")) +
  theme(
    plot.title = element_text(size = 8),        # Title font size
    axis.title = element_text(size = 6),        # Axis titles (UMAP_1, UMAP_2)
    axis.text = element_text(size = 6),         # Axis tick labels (numbers)
    legend.title = element_text(size = 6),      # Legend title font size
    legend.text = element_text(size = 5)        # Legend labels font size
  ) +
  labs(color = "Segmn Mean")                    # Custom legend title
dev.off()

# ✅ The segm.mean values in your SEG file represent the average (or consensus) copy number signal 
# for all cells assigned to a single clone.

# The distribution of per-cell mean CNA values (segm.mean)
# across the genomic interval:
# chr3:158,644,278–194,498,364
# only in cells that had a defined value (no NA)

# 🔬 Biological Interpretation:
# Higher values (e.g., > 0.2) = cells with copy number gain in that chr3 region
# Lower or negative values = cells with neutral or deleted signal
# Spread/shape of the violins = within-cluster heterogeneity in CNA signal


# ✅ The per-cell copy number signals — including segm.mean-like values — 
# are stored in the matrix file, typically called CNA_mtx_relat in SCEVAN.

# 🧬 To clarify:
# File/Object What it contains  Per Cell? Per Clone?
# CNA_mtx_relat Matrix of inferred copy number signals per segment/bin per cell (e.g., log2 ratios) ✅ Yes ❌ No
# SEG file (segm.mean)  Summary of copy number signal per segment aggregated across a clone ❌ No  ✅ Yes

######################################################################################
######################################################################################
######################################################################################
###################################################################################### display SegmMean profiles per chrom

# Use current directory for SEG files
seg_files <- list.files(".", pattern = "\\.seg$", full.names = TRUE)

# Initialize list to collect data
plot_list <- list()

# Loop through SEG files
for (file in seg_files) {
  clone_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the file
  seg_df <- read.table(file, header = TRUE)
  
  # Check expected columns
  if (all(c("Chr", "Pos", "End", "segm.mean") %in% colnames(seg_df))) {
    seg_df$clone <- clone_name
    seg_df$Midpoint <- (seg_df$Pos + seg_df$End) / 2
    plot_list[[clone_name]] <- seg_df
  } else {
    warning(paste("Skipping file due to missing columns:", file))
  }
}

# Combine all clones
seg_all <- bind_rows(plot_list)

# Filter for chromosome 3
# seg_chr3 <- seg_all %>% filter(Chr == 3)
seg_chr1 <- seg_all %>% filter(Chr == 1)


# Plot with enhancements
p <- ggplot(seg_chr1, aes(x = Midpoint / 1e6, y = segm.mean, color = clone)) +
  geom_line(size = 1.3) +  # Thicker lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.6) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "darkred", size = 0.6) +
  facet_wrap(~clone, ncol = 1, scales = "free_y") +
  labs(
    title = "Segm Mean profiles by clone",
    x = "Genomic Position (Mb)",
    y = "Segment Mean (log2R)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),      # Reduced from 12
    axis.title = element_text(size = 8),      # Reduced from 12
    plot.title = element_text(size = 8),      # Reduced from 14
    axis.text = element_text(size = 8),       # Added smaller axis text
    legend.text = element_text(size = 7),     # Legend text size
    legend.title = element_text(size = 8)     # Legend title size
  )

# Display plot
png("example_chr1_SegmMean_profiles.png", width = 1000, height = 1200, res = 150)
print(p)
dev.off()

######################################################################################
######################################################################################
######################################################################################
######################################################################################

# 🧬 Goal:
# To visualize clonal copy number alterations (CNAs) along chromosome 3 by reading multiple .seg files 
# — one per clone — and plotting their segm.mean (log2 ratio) values along the genomic coordinate.

######################################################################################
######################################################################################
######################################################################################
######################################################################################
###################################################################################### display CN profiles per chrom

# To show absolute copy number (CN) instead of segm.mean :

# Use current directory for SEG files
seg_files <- list.files(".", pattern = "\\.seg$", full.names = TRUE)

# Initialize list to collect data
plot_list <- list()

# Loop through SEG files
for (file in seg_files) {
  clone_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the file
  seg_df <- read.table(file, header = TRUE)
  
  # Check expected columns
  if (all(c("Chr", "Pos", "End", "CN") %in% colnames(seg_df))) {
    seg_df$clone <- clone_name
    seg_df$Midpoint <- (seg_df$Pos + seg_df$End) / 2
    plot_list[[clone_name]] <- seg_df
  } else {
    warning(paste("Skipping file due to missing columns:", file))
  }
}

# Combine all clones
seg_all <- bind_rows(plot_list)

# Filter for chromosome 3
# seg_chr3 <- seg_all %>% filter(Chr == 3)
seg_chr1 <- seg_all %>% filter(Chr == 1)

# Plot using CN instead of segm.mean
p2 <- ggplot(seg_chr1, aes(x = Midpoint / 1e6, y = CN, color = clone)) +
  geom_line(size = 1.3) +  # Thicker lines
  facet_wrap(~clone, ncol = 1, scales = "free_y") +
  labs(
    title = "Copy Number profiles by clone",
    x = "Genomic Position (Mb)",
    y = "Copy Number (CN)"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))

# Display plot
print(p2)

# Save to file
# ggsave("chr3_CN_by_clone.pdf", plot = p2, width = 10, height = 6 + length(unique(seg_chr3$clone)) * 1.5)

png("example_chr1_CNA_profiles.png", width = 1000, height = 1200, res = 150)
print(p2)
dev.off()

######################################################################################
######################################################################################
######################################################################################
######################################################################################

# Step 1: Read all SEG files

seg_files <- list.files(".", pattern = "\\.seg$", full.names = TRUE)
seg_list <- list()

for (file in seg_files) {
  clone_name <- tools::file_path_sans_ext(basename(file))
  seg <- read.table(file, header = TRUE)

  # Check required columns
  if (all(c("Chr", "Pos", "End", "CN", "segm.mean") %in% colnames(seg))) {
    seg$clone <- clone_name
    seg$Midpoint <- (seg$Pos + seg$End) / 2
    seg_list[[clone_name]] <- seg
  } else {
    warning(paste("Skipping file due to missing columns:", file))
  }
}

# Step 2: Combine and factor variables

all_seg <- bind_rows(seg_list)
all_seg$Chr <- factor(all_seg$Chr, levels = sort(unique(all_seg$Chr)))
all_seg$clone <- factor(all_seg$clone, levels = unique(all_seg$clone))

# Step 3: Plot CN profile

p_cn <- ggplot(all_seg, aes(x = Midpoint / 1e6, y = CN)) +
  geom_line(aes(group = 1), color = "steelblue", size = 0.8) +
  facet_grid(clone ~ Chr, scales = "free_x", space = "free_x") +
  labs(
    title = "Genome-wide Copy Number (CN) Profiles by Clone",
    x = "Genomic Position (Mb)",
    y = "Copy Number"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),           # Facet labels
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  # Changed from 6 to 8
    axis.text.y = element_text(size = 8),          # Y-axis text
    axis.title = element_text(size = 8),           # Axis titles
    plot.title = element_text(size = 8, face = "bold"),  # Changed from 14 to 8
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Step 4: Plot segm.mean profile
p_segmean <- ggplot(all_seg, aes(x = Midpoint / 1e6, y = segm.mean)) +
  geom_line(aes(group = 1), color = "firebrick", size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "darkred") +
  facet_grid(clone ~ Chr, scales = "free_x", space = "free_x") +
  labs(
    title = "Genome-wide Segment Mean (log2 ratio) Profiles by Clone",
    x = "Genomic Position (Mb)",
    y = "Segment Mean"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),           # Facet labels
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  # Changed from 6 to 8
    axis.text.y = element_text(size = 8),          # Y-axis text
    axis.title = element_text(size = 8),           # Axis titles
    plot.title = element_text(size = 8, face = "bold"),  # Changed from 14 to 8
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Step 5: Display both plots
print(p_cn)
print(p_segmean)

# Fixed the PNG outputs
png("genome_wide_segm_means_profiles_per_chromosome.png", width = 3200, height = 1200, res = 150)
print(p_segmean)
dev.off()

png("genome_wide_CN_profiles_per_chromosome.png", width = 3200, height = 1200, res = 150)
print(p_cn)
dev.off()

######################################################################################
###################################################################################### select a gene and display :
###################################################################################### Segm Mean
###################################################################################### Copy Number
###################################################################################### experimental code
# To visualize the segment mean (segm.mean) or Copy Number (CN) of a specific gene across all cells on a UMAP, 

# Set your gene of interest
# gene_of_interest <- "HES4"  # Change to the gene you re interested in

# Assumes the following are already loaded:
# - obj_scevan2: your Seurat object with UMAP
# - CNA_mtx_relat: matrix of segm.mean (rows = segments or genes, columns = cells)
# - count_mtx_annot: data frame with row-wise annotation for CNA_mtx_relat (must include gene name)

# Example structure of count_mtx_annot:
# head(count_mtx_annot)
# gene        chr     start        end
# TP53         17   7661779     7687550

# Set your gene of interest
gene_of_interest <- "HES4"  # Change to the gene you're interested in

# Assumes the following are already loaded:
# - obj_scevan2: your Seurat object with UMAP
# - CNA_mtx_relat: matrix of segm.mean (rows = segments or genes, columns = cells)
# - count_mtx_annot: data frame with row-wise annotation for CNA_mtx_relat (must include gene name)

# Step 1: Locate the row corresponding to your gene
gene_row <- which(count_mtx_annot$gene_name == gene_of_interest)
if (length(gene_row) == 0) {
  stop(paste("Gene", gene_of_interest, "not found in count_mtx_annot."))
}

# Step 2: Extract segm.mean or CN values for that gene across all cells
gene_cna <- CNA_mtx_relat[gene_row, ]  # returns a numeric vector

# Step 3: Ensure names match Seurat cells
# Step 4: Debug and fix the metadata assignment
seurat_cells <- Cells(obj_scevan2)
cna_cells <- names(gene_cna)

# Find intersection of cells
common_cells <- intersect(seurat_cells, cna_cells)
print(paste("Common cells:", length(common_cells)))

if (length(common_cells) == 0) {
  stop("No overlapping cells between Seurat object and CNA data")
}

# Subset Seurat object to only common cells
obj_subset <- subset(obj_scevan2, cells = common_cells)

# Get CNA values for common cells only
cna_values <- gene_cna[common_cells]
cna_values <- as.numeric(cna_values)
names(cna_values) <- NULL

# Add to subset metadata
obj_subset@meta.data[[paste0(gene_of_interest, "_CNA")]] <- cna_values

# Plot with color gradient from -0.5 to 0.5
p <- FeaturePlot(
  obj_subset,
  features = paste0(gene_of_interest, "_CNA"),
  pt.size = 1.2,  # Increased from 0.8 to 1.2
  order = TRUE
) + 
  scale_colour_gradient2(
    low = "#4575b4",      # Blue
    mid = "lightgray",    # Light gray
    high = "#d73027",     # Red
    midpoint = 0,
    breaks = c(-0.5, -0.2, 0, 0.2, 0.5),  # Updated breaks for -0.5 to 0.5 range
    labels = c("-0.5", "-0.2", "0.0", "0.2", "0.5"),  # Updated labels
    limits = c(-0.5, 0.5),  # Changed from c(-1, 1) to c(-0.5, 0.5)
    name = "Segm Mean"
  ) +
  ggtitle(paste(gene_of_interest)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 6),        # Reduced from 8 to 6
    legend.title = element_text(size = 6),       # Reduced from 8 to 6
    legend.key.size = unit(0.4, "cm"),          # Make color boxes smaller
    legend.key.height = unit(0.3, "cm"),        # Make legend height smaller
    legend.key.width = unit(0.3, "cm"),         # Make legend width smaller
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0), # Remove legend margins
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Save the plot
ggsave(
  filename = paste0("gene_", gene_of_interest, "_segm_mean_UMAP.png"),
  plot = p,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)

# print(p)
source("../script.plot.gene.segmn.mean.R")

######################################################################################
######################################################################################

# Load the function
source("../script.plot.gene.segmn.mean.function.R")

# Create plots for HES4
results <- create_gene_cna_plots(
  seurat_obj = obj_subset,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_name = "HES4",
  pt_size = 0.5,        # Small points
  figure_width = 6,     # Compact width
  figure_height = 5     # Compact height
)

unique(count_mtx_annot$gene_name)       # See all unique gene names
grep("TP53", count_mtx_annot$gene_name, value = TRUE, ignore.case = TRUE)  # Find matches ignoring case

# Create plots for multiple genes at once
gene_list <- c("HES4", "TP53", "MYC", "EGFR")
all_results <- create_multiple_gene_plots(
  seurat_obj = obj_subset,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_list = gene_list,
  pt_size = 0.5
)

# The co-deletion of 1p/19q is a defining feature of oligodendrogliomas and is associated with better prognosis 
# and responsiveness to therapy.

# While deletions on chromosome 6p are observed in various gliomas, specific genes affected on 6p in low-grade gliomas 
# are less well-characterized and may require further research.

# Genes deleted in low-grade gliomas (1p and 19q)

genes_1p <- c(
  "FUBP1", "NOTCH2", "CHD5", "PRKCG", "SLC17A7",
  "RIMS3", "KIAA1324", "AK5", "SLC6A17", "CD22",
  "HPCA", "MAG"
)

genes_19q <- c(
  "CIC", "NOTCH3", "CAPG", "PPP1R15A", "ZNF342",
  "NDUFA3", "RABAC1", "SLC25A5", "TUBB4A", "APOE",
  "BCL2L12", "RPS19"
)

# Genes amplified on chromosome 11q in gliomas
# 11q13 region : a locus known for amplification events in various cancers, including gliomas. 

genes_11q_amplified <- c("CCND1", "FGF3", "FGF4", "FGF19")

print(genes_1p)
print(genes_19q)
print(genes_11q_amplified)

all_results <- create_multiple_gene_plots(
  seurat_obj = obj_scevan2,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_list = genes_1p,
  pt_size = 0.5
)

all_results <- create_multiple_gene_plots(
  seurat_obj = obj_scevan2,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_list = genes_19q,
  pt_size = 0.5
)

all_results <- create_multiple_gene_plots(
  seurat_obj = obj_scevan2,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_list = genes_11q_amplified,
  pt_size = 0.5
)

######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################

# Load the categorical function
source("../script.plot.gene.copy.number.function.R")

# Create categorical plots for HES4
results <- create_categorical_gene_cna_plots(
  seurat_obj = obj_scevan2,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_name = "HES4",
  gain_threshold = 0.2,   # Customize gain threshold
  loss_threshold = -0.2,  # Customize loss threshold
  pt_size = 0.5
)

# Multiple genes
gene_list <- c("HES4", "TP53", "MYC")
all_results <- create_multiple_categorical_plots(
  seurat_obj = obj_scevan2,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_list = gene_list,
  gain_threshold = 0.2,  # More sensitive threshold
  loss_threshold = -0.2
)

# View category counts
# results$category_summary

all_results <- create_multiple_categorical_plots(
  seurat_obj = obj_scevan2,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_list = genes_1p,
  pt_size = 0.5
)

all_results <- create_multiple_categorical_plots(
  seurat_obj = obj_scevan2,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_list = genes_19q,
  pt_size = 0.5
)

all_results <- create_multiple_categorical_plots(
  seurat_obj = obj_scevan2,
  cna_matrix = CNA_mtx_relat,
  gene_annotations = count_mtx_annot,
  gene_list = genes_11q_amplified,
  pt_size = 0.5
)

######################################################################################
######################################################################################

save.image(file = paste0(file_base, ".RData"))

######################################################################################
######################################################################################
######################################################################################
