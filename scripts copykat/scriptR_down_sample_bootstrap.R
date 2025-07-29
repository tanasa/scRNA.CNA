#!/usr/bin/env Rscript

# Bootstrap Downsampling - 10 Iterations
# Automatically generates 10 downsampled datasets

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(knitr)
library(pander)
library(patchwork)
library(pheatmap)

# ===============================================
# LOAD DATA ONCE
# ===============================================

cat("=== LOADING DATA ===\n")
it = readRDS("merged_AG_afCB_afDF_res0.3_annotation_050625.rds")

# Basic dataset summary
cat("Dataset:", deparse(substitute(it)), "\n")
cat("Total cells:", ncol(it), "\n")
cat("Samples:", paste(unique(it@meta.data$orig.ident), collapse = ", "), "\n")
cat("Cell types:", paste(sort(unique(it@meta.data$celltype)), collapse = ", "), "\n")

# Show sample-wise totals
cat("\nSample totals:\n")
sample_counts <- table(it@meta.data$orig.ident)
print(sample_counts)

min_cells <- min(sample_counts)
cat("\nTarget size per sample:", min_cells, "cells\n")
cat("Total cells before:", sum(sample_counts), "\n")
cat("Total cells after per iteration:", min_cells * length(sample_counts), "\n\n")

# ===============================================
# ITERATE THROUGH 10 DOWNSAMPLING RUNS
# ===============================================

# Storage for summary statistics
all_results <- list()
summary_stats <- data.frame()

cat("=== STARTING 10 DOWNSAMPLING ITERATIONS ===\n\n")

for (i in 1:10) {
  
  cat("===============================================\n")
  cat("           DOWNSAMPLING RUN", i, "OF 10\n")
  cat("===============================================\n")
  
  # ===============================================
  # DOWNSAMPLE EACH SAMPLE
  # ===============================================
  
  cat("=== DOWNSAMPLING RUN", i, "===\n")
  
  # Set seed based on run number for different results each time
  set.seed(12345 + i)
  cat("Using seed:", 12345 + i, "\n")
  
  # Get selected cells from each sample
  selected_cells <- c()
  
  for (sample in names(sample_counts)) {
    sample_cells <- rownames(it@meta.data)[it@meta.data$orig.ident == sample]
    
    if (length(sample_cells) == min_cells) {
      selected <- sample_cells
      cat(sample, ": keeping all", length(selected), "cells\n")
    } else {
      selected <- sample(sample_cells, size = min_cells, replace = FALSE)
      cat(sample, ": selected", length(selected), "out of", length(sample_cells), "cells\n")
    }
    
    selected_cells <- c(selected_cells, selected)
  }
  
  cat("Total selected cells:", length(selected_cells), "\n")
  
  # ===============================================
  # CREATE DOWNSAMPLED OBJECT
  # ===============================================
  
  cat("\n=== CREATING DOWNSAMPLED OBJECT ===\n")
  
  it_downsampled <- subset(it, cells = selected_cells)
  
  cat("Original:", ncol(it), "cells\n")
  cat("Downsampled:", ncol(it_downsampled), "cells\n")
  
  # ===============================================
  # VERIFY RESULTS
  # ===============================================
  
  cat("\n=== VERIFICATION ===\n")
  new_counts <- table(it_downsampled@meta.data$orig.ident)
  print(new_counts)
  
  if (all(new_counts == min_cells)) {
    cat("✓ SUCCESS: All samples have", min_cells, "cells\n")
  } else {
    cat("✗ ERROR: Unequal sample sizes\n")
  }
  
  # ===============================================
  # SAVE RESULTS WITH DYNAMIC NAMING
  # ===============================================
  
  cat("\n=== SAVING RESULTS ===\n")
  
  # Create dynamic file names
  rds_filename <- paste0("down_sample_run_", i, ".rds")
  csv_filename <- paste0("down_sample_run_", i, "_selected_cells.csv")
  table_filename <- paste0("down_sample_run_", i, "_cell_distribution.csv")
  plot1_filename <- paste0("down_sample_run_", i, "_umap_celltypes.png")
  plot2_filename <- paste0("down_sample_run_", i, "_umap_samples.png")
  combined_filename <- paste0("down_sample_run_", i, "_umap_combined.png")
  png_filename <- paste0("down_sample_run_", i, "_heatmap.png")
  
  # Save downsampled object
  saveRDS(it_downsampled, rds_filename)
  cat("✓ Saved:", rds_filename, "\n")
  
  # Save selected cell list
  write.csv(data.frame(cell_barcode = selected_cells), 
            csv_filename, row.names = FALSE)
  cat("✓ Saved:", csv_filename, "\n")
  
  # ===============================================
  # CREATE SUMMARY TABLE
  # ===============================================
  
  # Create cross-tabulation table
  celltype_by_sample <- table(it_downsampled@meta.data$celltype, it_downsampled@meta.data$orig.ident)
  counts_df <- as.data.frame.matrix(celltype_by_sample)
  percentages_df <- round(prop.table(as.matrix(counts_df), margin = 2) * 100, 1)
  
  # Store percentages for stability analysis
  all_results[[paste0("run_", i)]] <- percentages_df
  
  # Combined format (Count (Percentage%))
  combined_df <- counts_df
  for(row in 1:nrow(counts_df)) {
    for(col in 1:ncol(counts_df)) {
      combined_df[row, col] <- paste0(counts_df[row, col], " (", percentages_df[row, col], "%)")
    }
  }
  
  # Save the table
  write.csv(combined_df, table_filename, row.names = TRUE)
  cat("✓ Saved table:", table_filename, "\n")
  
  # ===============================================
  # VISUALIZATION
  # ===============================================
  
  cat("\n=== CREATING VISUALIZATIONS ===\n")
  
  # UMAP plots if available
  if ("umap" %in% names(it_downsampled@reductions)) {
    cat("Creating UMAP plots...\n")
    
    # UMAP by cell type
    p1 <- DimPlot(it_downsampled, reduction = "umap", group.by = "celltype", 
                  label = TRUE, pt.size = 0.5) +
      labs(title = paste("UMAP: Cell Types - Run", i)) +
      theme(legend.position = "bottom")
    
    # UMAP by sample
    p2 <- DimPlot(it_downsampled, reduction = "umap", group.by = "orig.ident", 
                  pt.size = 0.5) +
      labs(title = paste("UMAP: Samples - Run", i)) +
      theme(legend.position = "bottom")
    
    # Display plots side by side (only for first run to avoid cluttering output)
    if (i == 1) {
      options(repr.plot.width = 16, repr.plot.height = 6)
      combined_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(6, 6))
      print(combined_plot)
    }
    
    # Save individual plots
    ggsave(plot1_filename, p1, width = 10, height = 8, dpi = 300)
    ggsave(plot2_filename, p2, width = 8, height = 8, dpi = 300)
    
    # Save combined plot
    combined_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(6, 6))
    ggsave(combined_filename, combined_plot, width = 20, height = 8, dpi = 300)
    
    cat("✓ Saved plots:", plot1_filename, "and", plot2_filename, "\n")
    cat("✓ Saved combined plot:", combined_filename, "\n")
  }
  
  # Heatmap
  heatmap_data <- as.matrix(percentages_df)
  
  png(png_filename, width = 10, height = 8, units = "in", res = 300)
  
  pheatmap(heatmap_data,
           main = paste("Cell Type Proportions by Sample (%) - Run", i),
           display_numbers = TRUE,
           number_format = "%.1f",
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("white", "red"))(50))
  
  dev.off()
  cat("✓ Saved heatmap:", png_filename, "\n")
  
  # ===============================================
  # CELL TYPE STATISTICS FOR THIS RUN
  # ===============================================
  
  # Calculate overall cell type percentages for this run
  overall_celltype <- table(it_downsampled@meta.data$celltype)
  overall_pct <- round(overall_celltype / sum(overall_celltype) * 100, 2)
  
  # Store in summary stats
  run_stats <- data.frame(
    Run = i,
    Seed = 12345 + i,
    Total_Cells = ncol(it_downsampled),
    CellTypes = paste(names(overall_pct), collapse = ";"),
    Percentages = paste(overall_pct, collapse = ";"),
    stringsAsFactors = FALSE
  )
  
  summary_stats <- rbind(summary_stats, run_stats)
  
  cat("\n=== RUN", i, "COMPLETE ===\n")
  cat("Files created for run", i, ":\n")
  cat("- Seurat object:", rds_filename, "\n")
  cat("- Selected cells:", csv_filename, "\n") 
  cat("- Distribution table:", table_filename, "\n")
  cat("- UMAP plots:", plot1_filename, ",", plot2_filename, "\n")
  cat("- Combined UMAP:", combined_filename, "\n")
  cat("- Heatmap:", png_filename, "\n\n")
  
}

# ===============================================
# STABILITY ANALYSIS ACROSS ALL RUNS
# ===============================================

cat("===============================================\n")
cat("           STABILITY ANALYSIS\n")
cat("===============================================\n")

# Combine all results for stability analysis
cat("=== ANALYZING STABILITY ACROSS 10 RUNS ===\n")

# Get all unique cell types
all_celltypes <- unique(unlist(lapply(all_results, rownames)))

# Create matrix to store all percentages
stability_matrix <- matrix(0, nrow = 10, ncol = length(all_celltypes))
colnames(stability_matrix) <- all_celltypes
rownames(stability_matrix) <- paste0("Run_", 1:10)

# Fill the matrix
for (run_idx in 1:10) {
  run_data <- all_results[[run_idx]]
  
  # For each cell type, get the mean percentage across all samples
  for (celltype in rownames(run_data)) {
    if (celltype %in% all_celltypes) {
      # Calculate mean percentage across samples for this cell type
      mean_pct <- mean(as.numeric(run_data[celltype, ]))
      stability_matrix[run_idx, celltype] <- mean_pct
    }
  }
}

# Calculate stability statistics
stability_stats <- data.frame(
  CellType = colnames(stability_matrix),
  Mean_Percent = round(colMeans(stability_matrix), 2),
  SD_Percent = round(apply(stability_matrix, 2, sd), 2),
  Min_Percent = round(apply(stability_matrix, 2, min), 2),
  Max_Percent = round(apply(stability_matrix, 2, max), 2),
  CV_Percent = round(apply(stability_matrix, 2, sd) / colMeans(stability_matrix) * 100, 1),
  stringsAsFactors = FALSE
)

# Sort by mean percentage
stability_stats <- stability_stats[order(stability_stats$Mean_Percent, decreasing = TRUE), ]

cat("Cell type stability across 10 iterations:\n")
print(kable(stability_stats, 
            caption = "Cell Type Stability Statistics Across 10 Runs",
            col.names = c("Cell Type", "Mean %", "SD", "Min %", "Max %", "CV %")))

# Identify stable vs variable cell types
stable_celltypes <- stability_stats[stability_stats$CV_Percent < 10 & stability_stats$Mean_Percent > 1, ]
variable_celltypes <- stability_stats[stability_stats$CV_Percent > 25 & stability_stats$Mean_Percent > 1, ]

cat("\nSTABLE cell types (CV < 10%, Mean > 1%):\n")
if (nrow(stable_celltypes) > 0) {
  for (i in 1:nrow(stable_celltypes)) {
    cat("  ", stable_celltypes$CellType[i], ": ", stable_celltypes$Mean_Percent[i], 
        "% ± ", stable_celltypes$SD_Percent[i], " (CV: ", stable_celltypes$CV_Percent[i], "%)\n")
  }
} else {
  cat("  No highly stable cell types found\n")
}

cat("\nVARIABLE cell types (CV > 25%, Mean > 1%):\n")
if (nrow(variable_celltypes) > 0) {
  for (i in 1:nrow(variable_celltypes)) {
    cat("  ", variable_celltypes$CellType[i], ": ", variable_celltypes$Mean_Percent[i], 
        "% ± ", variable_celltypes$SD_Percent[i], " (CV: ", variable_celltypes$CV_Percent[i], "%)\n")
  }
} else {
  cat("  No highly variable cell types found\n")
}

# ===============================================
# SAVE COMPREHENSIVE RESULTS
# ===============================================

cat("\n=== SAVING COMPREHENSIVE RESULTS ===\n")

# Save run summary
write.csv(summary_stats, "down_sampling_runs_summary.csv", row.names = FALSE)
cat("✓ Saved: down_sampling_runs_summary.csv\n")

# Save stability analysis
write.csv(stability_stats, "down_sampling_celltype_stability_analysis.csv", row.names = FALSE)
cat("✓ Saved: down_sampling_celltype_stability_analysis.csv\n")

# Save stability matrix
write.csv(stability_matrix, "down_sampling_celltype_percentages_all_runs.csv", row.names = TRUE)
cat("✓ Saved: down_sampling_celltype_percentages_all_runs.csv\n")

# ===============================================
# SAVE STABLE AND VARIABLE CELL TYPES
# ===============================================

cat("\n=== SAVING STABLE AND VARIABLE CELL TYPES ===\n")

# Save stable cell types
if (nrow(stable_celltypes) > 0) {
  stable_celltypes$Classification <- "STABLE (CV < 10%, Mean > 1%)"
  stable_celltypes$Interpretation <- "Low variability across runs - Reliable for analysis"
  
  write.csv(stable_celltypes, "down_sampling_STABLE_celltypes.csv", row.names = FALSE)
  cat("✓ Saved: down_sampling_STABLE_celltypes.csv\n")
  cat("  - Contains", nrow(stable_celltypes), "stable cell types\n")
} else {
  cat("! No stable cell types found (CV < 10%, Mean > 1%)\n")
}

# Save variable cell types
if (nrow(variable_celltypes) > 0) {
  variable_celltypes$Classification <- "VARIABLE (CV > 25%, Mean > 1%)"
  variable_celltypes$Interpretation <- "High variability across runs - Use with caution"
  
  write.csv(variable_celltypes, "down_sampling_VARIABLE_celltypes.csv", row.names = FALSE)
  cat("✓ Saved: down_sampling_VARIABLE_celltypes.csv\n")
  cat("  - Contains", nrow(variable_celltypes), "variable cell types\n")
} else {
  cat("! No highly variable cell types found (CV > 25%, Mean > 1%)\n")
}

# Save combined classification file
all_classified <- rbind(
  if(nrow(stable_celltypes) > 0) stable_celltypes else NULL,
  if(nrow(variable_celltypes) > 0) variable_celltypes else NULL
)

if (!is.null(all_classified) && nrow(all_classified) > 0) {
  write.csv(all_classified, "down_sampling_classified_celltypes.csv", row.names = FALSE)
  cat("✓ Saved: down_sampling_classified_celltypes.csv\n")
  cat("  - Contains both stable and variable cell types with classifications\n")
}

# Create a summary classification report
cat("\n=== CREATING CLASSIFICATION SUMMARY REPORT ===\n")

classification_summary <- c(
  "BOOTSTRAP DOWNSAMPLING CLASSIFICATION REPORT",
  paste("Generated on:", Sys.time()),
  "",
  "CLASSIFICATION CRITERIA:",
  "- STABLE: CV < 10% AND Mean > 1%",
  "  → Low variability across runs, reliable for analysis",
  "- VARIABLE: CV > 25% AND Mean > 1%", 
  "  → High variability across runs, use with caution",
  "- MODERATE: 10% ≤ CV ≤ 25% AND Mean > 1%",
  "  → Moderate variability, acceptable for analysis",
  "",
  "RESULTS:",
  paste("- Total cell types analyzed:", length(all_celltypes)),
  paste("- STABLE cell types:", nrow(stable_celltypes)),
  paste("- VARIABLE cell types:", nrow(variable_celltypes)),
  paste("- MODERATE cell types:", sum(stability_stats$CV_Percent >= 10 & stability_stats$CV_Percent <= 25 & stability_stats$Mean_Percent > 1)),
  paste("- Rare cell types (<1% mean):", sum(stability_stats$Mean_Percent <= 1)),
  "",
  "STABLE CELL TYPES (Recommended for analysis):"
)

if (nrow(stable_celltypes) > 0) {
  for (i in 1:nrow(stable_celltypes)) {
    line <- paste("- ", stable_celltypes$CellType[i], ": ", 
                  stable_celltypes$Mean_Percent[i], "% ± ", 
                  stable_celltypes$SD_Percent[i], " (CV: ", 
                  stable_celltypes$CV_Percent[i], "%)")
    classification_summary <- c(classification_summary, line)
  }
} else {
  classification_summary <- c(classification_summary, "- None found")
}

classification_summary <- c(classification_summary, "", "VARIABLE CELL TYPES (Use with caution):")

if (nrow(variable_celltypes) > 0) {
  for (i in 1:nrow(variable_celltypes)) {
    line <- paste("- ", variable_celltypes$CellType[i], ": ", 
                  variable_celltypes$Mean_Percent[i], "% ± ", 
                  variable_celltypes$SD_Percent[i], " (CV: ", 
                  variable_celltypes$CV_Percent[i], "%)")
    classification_summary <- c(classification_summary, line)
  }
} else {
  classification_summary <- c(classification_summary, "- None found")
}

classification_summary <- c(classification_summary, 
  "",
  "RECOMMENDATIONS:",
  "1. Use STABLE cell types for differential expression analysis",
  "2. Exercise caution with VARIABLE cell types",
  "3. Consider filtering out rare cell types (<1% mean abundance)",
  "4. Validate findings using multiple downsampling iterations",
  "",
  "FILES GENERATED:",
  "- down_sampling_STABLE_celltypes.csv",
  "- down_sampling_VARIABLE_celltypes.csv", 
  "- down_sampling_classified_celltypes.csv",
  "- down_sampling_classification_report.txt (this file)"
)

# Save the report
writeLines(classification_summary, "down_sampling_classification_report.txt")
cat("✓ Saved: down_sampling_classification_report.txt\n")
cat("  - Contains detailed classification summary and recommendations\n")

# ===============================================
# FINAL SUMMARY
# ===============================================

cat("\n===============================================\n")
cat("           FINAL SUMMARY\n")
cat("===============================================\n")

cat("✓ Completed 10 downsampling iterations successfully\n")
cat("✓ Each iteration:", min_cells, "cells per sample,", length(sample_counts), "samples\n")
cat("✓ Total cells per iteration:", min_cells * length(sample_counts), "\n")
cat("✓ Seeds used: 12346 to 12355\n")
cat("✓ Cell types analyzed:", length(all_celltypes), "\n")
cat("✓ Stable cell types (CV < 10%):", nrow(stable_celltypes), "\n")
cat("✓ Variable cell types (CV > 25%):", nrow(variable_celltypes), "\n")

cat("\nFiles generated:\n")
cat("- 10 Seurat objects: down_sample_run_1.rds to down_sample_run_10.rds\n")
cat("- 10 Cell lists: down_sample_run_X_selected_cells.csv\n")
cat("- 10 Distribution tables: down_sample_run_X_cell_distribution.csv\n")
cat("- 30 UMAP plots: Individual and combined for each run\n")
cat("- 10 Heatmaps: down_sample_run_X_heatmap.png\n")
cat("- 3 Summary files: downsampling_runs_summary.csv, celltype_stability_analysis.csv, celltype_percentages_all_runs.csv\n")

cat("\nRecommendation: Use celltype_stability_analysis.csv to identify\n")
cat("robust cell type differences that are consistent across multiple downsampling iterations.\n")

cat("\n=== ALL 10 ITERATIONS COMPLETE ===\n")
