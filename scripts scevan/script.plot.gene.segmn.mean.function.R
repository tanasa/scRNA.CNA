# comprehensive_gene_cna_plots.R
# Complete function to create all CNA plots for a gene

# Suppress ComplexHeatmap messages
library(ComplexHeatmap) 
ht_opt$message = FALSE

# ===================================================

create_gene_cna_plots <- function(seurat_obj, cna_matrix, gene_annotations, gene_name, 
                                 pt_size = 0.5, figure_width = 6, figure_height = 5) {
  
  cat("Creating comprehensive CNA plots for", gene_name, "...\n")
  
  # ============================================================================
  # Step 1: Data preparation and validation with improved gene searching
  # ============================================================================
  
  # Check if gene exists - exact match first
  gene_row <- which(gene_annotations$gene_name == gene_name)
  
  if (length(gene_row) == 0) {
    # Try case-insensitive search
    cat("Exact match not found. Trying case-insensitive search...\n")
    gene_matches <- grep(paste0("^", gene_name, "$"), gene_annotations$gene_name, 
                        value = TRUE, ignore.case = TRUE)
    
    if (length(gene_matches) > 0) {
      cat("Found case-insensitive match:", gene_matches[1], "\n")
      gene_row <- which(gene_annotations$gene_name == gene_matches[1])
      gene_name <- gene_matches[1]  # Update to correct case
    } else {
      # Try partial match
      cat("No exact match found. Searching for partial matches...\n")
      partial_matches <- grep(gene_name, gene_annotations$gene_name, 
                             value = TRUE, ignore.case = TRUE)
      
      if (length(partial_matches) > 0) {
        cat("Found", length(partial_matches), "partial matches:\n")
        print(head(partial_matches, 10))
        cat("Please use exact gene name. Suggested:", partial_matches[1], "\n")
        stop(paste("Gene", gene_name, "not found. Try:", partial_matches[1]))
      } else {
        # No matches found - provide helpful info
        cat("No matches found for", gene_name, "\n")
        cat("Total genes available:", length(unique(gene_annotations$gene_name)), "\n")
        cat("First 10 gene names as examples:\n")
        print(head(unique(gene_annotations$gene_name), 10))
        stop(paste("Gene", gene_name, "not found in gene annotations."))
      }
    }
  }
  
  cat("Found gene", gene_name, "at row", gene_row, "\n")
  
  # Extract CNA values for the gene
  gene_cna <- cna_matrix[gene_row, ]
  
  # Find common cells between Seurat object and CNA data
  seurat_cells <- Cells(seurat_obj)
  cna_cells <- names(gene_cna)
  common_cells <- intersect(seurat_cells, cna_cells)
  
  cat("Common cells found:", length(common_cells), "\n")
  
  if (length(common_cells) == 0) {
    stop("No overlapping cells between Seurat object and CNA data")
  }
  
  # Subset Seurat object to only common cells
  obj_subset <- subset(seurat_obj, cells = common_cells)
  
  # Get CNA values for common cells only
  cna_values <- gene_cna[common_cells]
  cna_values <- as.numeric(cna_values)
  names(cna_values) <- NULL
  
  # Add to subset metadata
  obj_subset@meta.data[[paste0(gene_name, "_CNA")]] <- cna_values
  
  # Also add to original Seurat object for compatibility
  seurat_obj <- AddMetaData(seurat_obj, 
                           metadata = setNames(list(cna_values), paste0(gene_name, "_CNA")),
                           col.name = paste0(gene_name, "_CNA"))
  
  # Print data summary
  cat("Data summary for", gene_name, ":\n")
  cat("Min:", round(min(cna_values, na.rm = TRUE), 3), "\n")
  cat("Max:", round(max(cna_values, na.rm = TRUE), 3), "\n")
  cat("Mean:", round(mean(cna_values, na.rm = TRUE), 3), "\n")
  cat("Median:", round(median(cna_values, na.rm = TRUE), 3), "\n")
  
  # ============================================================================
  # Plot 1: Custom gradient with -0.5 to 0.5 range
  # ============================================================================
  
  p1 <- FeaturePlot(
    obj_subset,
    features = paste0(gene_name, "_CNA"),
    pt.size = pt_size,
    order = TRUE
  ) + 
    scale_colour_gradient2(
      low = "#4575b4",      # Blue
      mid = "lightgray",    # Light gray
      high = "#d73027",     # Red
      midpoint = 0,
      breaks = c(-0.5, -0.2, 0, 0.2, 0.5),
      labels = c("-0.5", "-0.2", "0.0", "0.2", "0.5"),
      limits = c(-0.5, 0.5),
      name = "Segm Mean"
    ) +
    ggtitle(paste(gene_name)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  print(p1)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_segm_mean_custom_range.png"),
    plot = p1,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 2: Blue-Red (RdBu palette)
  # ============================================================================
  
  p2 <- FeaturePlot(
    obj_subset,
    features = paste0(gene_name, "_CNA"),
    pt.size = pt_size
  ) + 
    scale_colour_distiller(
      palette = "RdBu", 
      direction = -1,
      name = "Segm Mean"
    ) +
    ggtitle(paste(gene_name)) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p2)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_segm_mean_RdBu.png"),
    plot = p2,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 3: Red-Yellow-Blue (RdYlBu palette)
  # ============================================================================
  
  p3 <- FeaturePlot(
    obj_subset,
    features = paste0(gene_name, "_CNA"),
    pt.size = pt_size
  ) + 
    scale_colour_distiller(
      palette = "RdYlBu", 
      direction = -1,
      name = "Segm Mean"
    ) +
    ggtitle(paste(gene_name)) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p3)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_segm_mean_red_yellow_blue.png"),
    plot = p3,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 4: Spectral palette
  # ============================================================================
  
  p4 <- FeaturePlot(
    obj_subset,
    features = paste0(gene_name, "_CNA"),
    pt.size = pt_size
  ) + 
    scale_colour_distiller(
      palette = "Spectral", 
      direction = -1,
      name = "Segm Mean"
    ) +
    ggtitle(paste(gene_name))+
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p4)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_segm_mean_spectral.png"),
    plot = p4,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 5: Custom gradient (autoscale)
  # ============================================================================
  
  p5 <- FeaturePlot(
    obj_subset,
    features = paste0(gene_name, "_CNA"),
    pt.size = pt_size
  ) + 
    scale_colour_gradientn(
      colors = c("darkblue", "blue", "lightblue", "white", "lightcoral", "red", "darkred"),
      name = "Segm Mean"
    ) +
    ggtitle(paste(gene_name)) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p5)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_segm_mean_custom_autoscale.png"),
    plot = p5,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 6: Custom gradient with fixed scale (-1 to +1)
  # ============================================================================
  
  p6 <- FeaturePlot(
    obj_subset,
    features = paste0(gene_name, "_CNA"),
    pt.size = pt_size
  ) + 
    scale_colour_gradientn(
      colors = c("darkblue", "blue", "lightblue", "white", "lightcoral", "red", "darkred"),
      values = scales::rescale(c(-1, -0.5, -0.2, 0, 0.2, 0.5, 1)),
      limits = c(-1, 1),
      name = "Segm Mean"
    ) +
    ggtitle(paste(gene_name)) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p6)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_segm_mean_custom_fixed_scale.png"),
    plot = p6,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Summary and return
  # ============================================================================
  
  cat("\nCompleted! Created 6 plots for", gene_name, ":\n")
  cat("1. gene_", gene_name, "_segm_mean_custom_range.png\n")
  cat("2. gene_", gene_name, "_segm_mean_blue_red.png\n")
  cat("3. gene_", gene_name, "_segm_mean_red_yellow_blue.png\n") 
  cat("4. gene_", gene_name, "_segm_mean_spectral.png\n")
  cat("5. gene_", gene_name, "_segm_mean_custom_autoscale.png\n")
  cat("6. gene_", gene_name, "_segm_mean_custom_fixed_scale.png\n")
  
  # Return a list with all plots and the updated objects
  return(list(
    plots = list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6),
    seurat_obj = seurat_obj,
    obj_subset = obj_subset,
    cna_values = cna_values,
    data_summary = list(
      min = min(cna_values, na.rm = TRUE),
      max = max(cna_values, na.rm = TRUE),
      mean = mean(cna_values, na.rm = TRUE),
      median = median(cna_values, na.rm = TRUE)
    )
  ))
}

# ============================================================================
# Helper functions for gene searching
# ============================================================================

# Function to search for genes with suggestions
search_gene <- function(gene_annotations, search_term) {
  cat("Searching for:", search_term, "\n")
  
  # Exact match
  exact_matches <- gene_annotations$gene_name[gene_annotations$gene_name == search_term]
  if (length(exact_matches) > 0) {
    cat("Exact match found:", exact_matches[1], "\n")
    return(exact_matches[1])
  }
  
  # Case-insensitive exact match
  case_matches <- grep(paste0("^", search_term, "$"), gene_annotations$gene_name, 
                      value = TRUE, ignore.case = TRUE)
  if (length(case_matches) > 0) {
    cat("Case-insensitive match found:", case_matches[1], "\n")
    return(case_matches[1])
  }
  
  # Partial matches
  partial_matches <- grep(search_term, gene_annotations$gene_name, 
                         value = TRUE, ignore.case = TRUE)
  if (length(partial_matches) > 0) {
    cat("Found", length(partial_matches), "partial matches:\n")
    print(head(partial_matches, 10))
    return(NULL)
  }
  
  cat("No matches found for", search_term, "\n")
  return(NULL)
}

# Function to list available genes
list_available_genes <- function(gene_annotations, pattern = NULL) {
  all_genes <- unique(gene_annotations$gene_name)
  
  if (is.null(pattern)) {
    cat("Total genes available:", length(all_genes), "\n")
    cat("First 20 genes:\n")
    print(head(all_genes, 20))
  } else {
    matches <- grep(pattern, all_genes, value = TRUE, ignore.case = TRUE)
    cat("Genes matching pattern '", pattern, "':\n")
    print(matches)
  }
  
  return(invisible(all_genes))
}

# Function to create plots for multiple genes
create_multiple_gene_plots <- function(seurat_obj, cna_matrix, gene_annotations, 
                                      gene_list, pt_size = 0.5, 
                                      figure_width = 6, figure_height = 5) {
  
  results <- list()
  
  for (gene in gene_list) {
    cat("\n", paste(rep("=", 50), collapse = ""), "\n")
    cat("Processing gene:", gene, "\n")
    cat(paste(rep("=", 50), collapse = ""), "\n")
    
    tryCatch({
      result <- create_gene_cna_plots(
        seurat_obj = seurat_obj,
        cna_matrix = cna_matrix,
        gene_annotations = gene_annotations,
        gene_name = gene,
        pt_size = pt_size,
        figure_width = figure_width,
        figure_height = figure_height
      )
      results[[gene]] <- result
      seurat_obj <- result$seurat_obj  # Update with new metadata
    }, error = function(e) {
      cat("Error processing", gene, ":", e$message, "\n")
    })
  }
  
  return(results)
}

# ============================================================================
# Usage instructions (commented out)
# ============================================================================

# How to use this function:
# 
# # Load the function
# source("comprehensive_gene_cna_plots.R")
# 
# # Search for genes before plotting
# search_gene(count_mtx_annot, "C11")           # Find C11 genes
# search_gene(count_mtx_annot, "orf30")         # Find orf30 genes
# list_available_genes(count_mtx_annot, "C11")  # List all C11 genes
# 
# # Create plots for a single gene (HES4)
# results <- create_gene_cna_plots(
#   seurat_obj = obj_scevan2,
#   cna_matrix = CNA_mtx_relat,
#   gene_annotations = count_mtx_annot,
#   gene_name = "HES4",
#   pt_size = 0.5,        # Point size
#   figure_width = 6,     # Figure width
#   figure_height = 5     # Figure height
# )
# 
# # Create plots for multiple genes
# gene_list <- c("HES4", "TP53", "MYC", "EGFR")
# all_results <- create_multiple_gene_plots(
#   seurat_obj = obj_scevan2,
#   cna_matrix = CNA_mtx_relat,
#   gene_annotations = count_mtx_annot,
#   gene_list = gene_list,
#   pt_size = 0.5,
#   figure_width = 6,
#   figure_height = 5
# )
# 
# # Access individual plots
# results$plots$p1  # Custom range plot
# results$plots$p2  # Blue-red plot
# # etc.
# 
# # Access data summary
# results$data_summary

cat("Comprehensive gene CNA plotting function loaded successfully!\n")
cat("Available functions:\n")
cat("- create_gene_cna_plots(): Create plots for single gene\n")
cat("- create_multiple_gene_plots(): Create plots for multiple genes\n") 
cat("- search_gene(): Search for gene names with suggestions\n")
cat("- list_available_genes(): List available genes with optional pattern\n")