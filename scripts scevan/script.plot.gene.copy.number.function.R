# categorical_gene_cna_plots.R
# Function to create categorical CNA plots (Gain/Loss/Neutral)

# Suppress ComplexHeatmap messages
library(ComplexHeatmap) 
ht_opt$message = FALSE

# ===========================================================

create_categorical_gene_cna_plots <- function(seurat_obj, cna_matrix, gene_annotations, gene_name, 
                                             gain_threshold = 0.2, loss_threshold = -0.2,
                                             pt_size = 0.5, figure_width = 6, figure_height = 5) {
  
  cat("Creating categorical CNA plots for", gene_name, "...\n")
  cat("Thresholds: Gain >", gain_threshold, ", Loss <", loss_threshold, "\n")
  
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
  
  # Create categorical classification
  library(dplyr)
  cna_categories <- case_when(
    cna_values > gain_threshold ~ "Gain",
    cna_values < loss_threshold ~ "Loss",
    TRUE ~ "Neutral"
  )
  
  # Add both continuous and categorical data to metadata
  obj_subset@meta.data[[paste0(gene_name, "_CNA")]] <- cna_values
  obj_subset@meta.data[[paste0(gene_name, "_CNA_class")]] <- cna_categories
  
  # Print data summary
  cat("Data summary for", gene_name, ":\n")
  cat("Min:", round(min(cna_values, na.rm = TRUE), 3), "\n")
  cat("Max:", round(max(cna_values, na.rm = TRUE), 3), "\n")
  cat("Mean:", round(mean(cna_values, na.rm = TRUE), 3), "\n")
  cat("Median:", round(median(cna_values, na.rm = TRUE), 3), "\n")
  
  # Print categorical summary
  cat("\nCategorical summary:\n")
  print(table(cna_categories))
  
  # ============================================================================
  # Plot 1: Categorical plot with custom colors
  # ============================================================================
  
  # Define colors for categories
  category_colors <- c("Loss" = "#4575b4", "Neutral" = "lightgray", "Gain" = "#d73027")
  
  p1 <- DimPlot(
    obj_subset,
    group.by = paste0(gene_name, "_CNA_class"),
    cols = category_colors,
    pt.size = pt_size,
    order = c("Gain", "Loss")  # Plot Gain and Loss on top
  ) + 
    ggtitle(paste(gene_name, "CNA Categories")) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p1)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_CNA_categorical.png"),
    plot = p1,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 2: Categorical plot with different color scheme
  # ============================================================================
  
  category_colors2 <- c("Loss" = "darkblue", "Neutral" = "gray90", "Gain" = "darkred")
  
  p2 <- DimPlot(
    obj_subset,
    group.by = paste0(gene_name, "_CNA_class"),
    cols = category_colors2,
    pt.size = pt_size,
    order = c("Gain", "Loss")
  ) + 
    ggtitle(paste(gene_name, "CNA (Dark Colors)")) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p2)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_CNA_categorical_dark.png"),
    plot = p2,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 3: Show only Gains
  # ============================================================================
  
  gain_colors <- c("Loss" = "lightgray", "Neutral" = "lightgray", "Gain" = "#d73027")
  
  p3 <- DimPlot(
    obj_subset,
    group.by = paste0(gene_name, "_CNA_class"),
    cols = gain_colors,
    pt.size = pt_size,
    order = c("Gain")
  ) + 
    ggtitle(paste(gene_name, "Gains Highlighted")) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p3)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_CNA_gains_only.png"),
    plot = p3,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 4: Show only Losses
  # ============================================================================
  
  loss_colors <- c("Loss" = "#4575b4", "Neutral" = "lightgray", "Gain" = "lightgray")
  
  p4 <- DimPlot(
    obj_subset,
    group.by = paste0(gene_name, "_CNA_class"),
    cols = loss_colors,
    pt.size = pt_size,
    order = c("Loss")
  ) + 
    ggtitle(paste(gene_name, "Losses Highlighted")) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p4)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_CNA_losses_only.png"),
    plot = p4,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 5: Continuous values for comparison
  # ============================================================================
  
  p5 <- FeaturePlot(
    obj_subset,
    features = paste0(gene_name, "_CNA"),
    pt.size = pt_size
  ) + 
    scale_colour_gradient2(
      low = "#4575b4", 
      mid = "lightgray", 
      high = "#d73027",
      midpoint = 0,
      name = "Segm Mean"
    ) +
    ggtitle(paste(gene_name, "Continuous Values")) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  print(p5)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_CNA_continuous.png"),
    plot = p5,
    width = figure_width, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Plot 6: Side-by-side comparison
  # ============================================================================
  
  library(patchwork)
  
  p6 <- p1 + p5 + plot_layout(ncol = 2) +
    plot_annotation(title = paste(gene_name, "- Categorical vs Continuous"))
  
  print(p6)
  
  ggsave(
    filename = paste0("gene_", gene_name, "_CNA_comparison.png"),
    plot = p6,
    width = figure_width * 2, height = figure_height, dpi = 300, bg = "white"
  )
  
  # ============================================================================
  # Summary and return
  # ============================================================================
  
  cat("\nCompleted! Created 6 plots for", gene_name, ":\n")
  cat("1. gene_", gene_name, "_CNA_categorical.png\n")
  cat("2. gene_", gene_name, "_CNA_categorical_dark.png\n")
  cat("3. gene_", gene_name, "_CNA_gains_only.png\n") 
  cat("4. gene_", gene_name, "_CNA_losses_only.png\n")
  cat("5. gene_", gene_name, "_CNA_continuous.png\n")
  cat("6. gene_", gene_name, "_CNA_comparison.png\n")
  
  # Return a list with all plots and data
  return(list(
    plots = list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6),
    obj_subset = obj_subset,
    cna_values = cna_values,
    cna_categories = cna_categories,
    category_summary = table(cna_categories),
    data_summary = list(
      min = min(cna_values, na.rm = TRUE),
      max = max(cna_values, na.rm = TRUE),
      mean = mean(cna_values, na.rm = TRUE),
      median = median(cna_values, na.rm = TRUE),
      gain_threshold = gain_threshold,
      loss_threshold = loss_threshold
    )
  ))
}

# ============================================================================
# Helper functions for gene searching (same as before)
# ============================================================================

search_gene <- function(gene_annotations, search_term) {
  cat("Searching for:", search_term, "\n")
  
  exact_matches <- gene_annotations$gene_name[gene_annotations$gene_name == search_term]
  if (length(exact_matches) > 0) {
    cat("Exact match found:", exact_matches[1], "\n")
    return(exact_matches[1])
  }
  
  case_matches <- grep(paste0("^", search_term, "$"), gene_annotations$gene_name, 
                      value = TRUE, ignore.case = TRUE)
  if (length(case_matches) > 0) {
    cat("Case-insensitive match found:", case_matches[1], "\n")
    return(case_matches[1])
  }
  
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

# Function to create plots for multiple genes
create_multiple_categorical_plots <- function(seurat_obj, cna_matrix, gene_annotations, 
                                            gene_list, gain_threshold = 0.2, loss_threshold = -0.2,
                                            pt_size = 0.5, figure_width = 6, figure_height = 5) {
  
  results <- list()
  
  for (gene in gene_list) {
    cat("\n", paste(rep("=", 50), collapse = ""), "\n")
    cat("Processing gene:", gene, "\n")
    cat(paste(rep("=", 50), collapse = ""), "\n")
    
    tryCatch({
      result <- create_categorical_gene_cna_plots(
        seurat_obj = seurat_obj,
        cna_matrix = cna_matrix,
        gene_annotations = gene_annotations,
        gene_name = gene,
        gain_threshold = gain_threshold,
        loss_threshold = loss_threshold,
        pt_size = pt_size,
        figure_width = figure_width,
        figure_height = figure_height
      )
      results[[gene]] <- result
    }, error = function(e) {
      cat("Error processing", gene, ":", e$message, "\n")
      cat("Use search_gene(gene_annotations, '", gene, "') to find alternatives\n")
    })
  }
  
  return(results)
}

# ============================================================================
# Usage instructions
# ============================================================================

cat("Categorical gene CNA plotting function loaded successfully!\n")
cat("Available functions:\n")
cat("- create_categorical_gene_cna_plots(): Create categorical plots for single gene\n")
cat("- create_multiple_categorical_plots(): Create categorical plots for multiple genes\n") 
cat("- search_gene(): Search for gene names with suggestions\n")

# Example usage:

# results <- create_categorical_gene_cna_plots(
#   seurat_obj = obj_scevan2,
#   cna_matrix = CNA_mtx_relat,
#   gene_annotations = count_mtx_annot,
#   gene_name = "HES4",
#   gain_threshold = 0.2,   # Threshold for gains
#   loss_threshold = -0.2,  # Threshold for losses
#   pt_size = 0.5
# )
