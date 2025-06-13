# https://www.youtube.com/watch?v=dj3yMpGfRtI

# Clean CNA Heatmap Script with Multiple Distance Metrics
# Modified to use CNX_mtx_clone_annot instead of CNA_mtx_relat_annot
# Load required libraries
library(ComplexHeatmap)
ht_opt$message = FALSE

library(circlize)
library(RColorBrewer)
library(grid)

# ----------------------------
# Step 1: Prepare CNA matrix
# ----------------------------

# Extract genomic annotation columns (first 5 columns)
genomic_info <- CNA_mtx_clone_annot[, 1:5]

# Extract CNA values (columns 6 onwards - the cell data)
cna_matrix <- CNA_mtx_clone_annot[, 6:ncol(CNA_mtx_clone_annot)]

# Ensure we have matching cells between CNA matrix and metadata
common_cells <- intersect(colnames(cna_matrix), rownames(metadata))
cna_matrix <- cna_matrix[, common_cells]
metadata_filtered <- metadata[common_cells, ]

print(paste("Number of cells:", ncol(cna_matrix)))
print(paste("Number of genes/segments:", nrow(cna_matrix)))

# ----------------------------
# Step 2: Prepare chromosome separation and colors
# ----------------------------

# Create chromosome factor with proper ordering
chr_levels <- c(as.character(1:22), "X", "Y")
chromosomes <- factor(genomic_info$seqnames, levels = chr_levels)

# Check what chromosomes are actually present in your data
unique_chrs <- unique(as.character(chromosomes))
print(paste("Chromosomes in data:", paste(unique_chrs, collapse = ", ")))

# Create chromosome boundaries for boxes
chr_boundaries <- list()
chr_colors <- list()
chr_start <- 1

for(i in 1:length(unique_chrs)) {
  chr_name <- unique_chrs[i]
  chr_segments <- which(as.character(chromosomes) == chr_name)
  chr_end <- chr_start + length(chr_segments) - 1
  
  chr_boundaries[[chr_name]] <- c(chr_start, chr_end)
  chr_colors[[chr_name]] <- ifelse(i %% 2 == 1, "darkblue", "darkred")
  
  chr_start <- chr_end + 1
}

print("Chromosome boundaries:")
print(chr_boundaries)

# ----------------------------
# Step 3: Prepare row annotation (subclones only)
# ----------------------------

# Extract subclone information
cell_subclones <- metadata_filtered$subclone
names(cell_subclones) <- rownames(metadata_filtered)

# Check what subclones are present
unique_subclones <- unique(cell_subclones[!is.na(cell_subclones)])
print(paste("Subclones in data:", paste(unique_subclones, collapse = ", ")))

# Create color palette for subclones
if(length(unique_subclones) > 0) {
  if(length(unique_subclones) <= 8) {
    subclone_colors <- RColorBrewer::brewer.pal(max(3, length(unique_subclones)), "Set1")[1:length(unique_subclones)]
  } else {
    subclone_colors <- rainbow(length(unique_subclones))
  }
  names(subclone_colors) <- unique_subclones
  
  # Handle NA values for subclones
  if(any(is.na(cell_subclones))) {
    subclone_colors <- c(subclone_colors, "white")
    names(subclone_colors)[length(subclone_colors)] <- "NA"
    cell_subclones[is.na(cell_subclones)] <- "NA"
  }
} else {
  # If no subclones, create a simple color scheme
  subclone_colors <- c("white")
  names(subclone_colors) <- "NA"
  cell_subclones[is.na(cell_subclones)] <- "NA"
}

print("Subclone colors:")
print(subclone_colors)

# Create row annotation with subclone only
row_anno <- rowAnnotation(
  Subclone = cell_subclones,
  col = list(Subclone = subclone_colors),
  annotation_name_side = "top",
  border = TRUE,
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  annotation_width = unit(8, "mm")
)

# ----------------------------
# Step 4: Prepare color scale for CNA values
# ----------------------------

# Calculate data range
cna_values <- as.matrix(cna_matrix)
data_range <- range(cna_values, na.rm = TRUE)
print(paste("Data range:", round(data_range[1], 3), "to", round(data_range[2], 3)))

# Fixed color scale between -1 to +1
col_breaks <- c(-1, -0.5, -0.2, 0, 0.2, 0.5, 1)

col_fun <- colorRamp2(col_breaks, 
                      colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(n = length(col_breaks)))

print("Fixed color scale: -1 to +1")
print("Copy number thresholds:")
print("Deep deletion: -1")
print("Deletion: -0.5") 
print("Shallow deletion: -0.2")
print("Normal: 0")
print("Shallow amplification: +0.2")
print("Amplification: +0.5")
print("High amplification: +1")

# ----------------------------
# Step 5: Define multiple distance metrics clustering functions
# ----------------------------

print("Setting up multiple distance metrics...")

# Function to compute cosine distance
cosine_distance <- function(x, y) {
  1 - sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

# Advanced clustering function with multiple distance metrics
create_clustering_function <- function(distance_method = "euclidean") {
  function(x) {
    tryCatch({
      # Different distance calculations
      if (distance_method == "cosine") {
        # Custom cosine distance implementation
        n <- nrow(x)
        d_matrix <- matrix(0, n, n)
        for(i in 1:(n-1)) {
          for(j in (i+1):n) {
            d_matrix[i,j] <- d_matrix[j,i] <- cosine_distance(x[i,], x[j,])
          }
        }
        d <- as.dist(d_matrix)
      } else if (distance_method == "minkowski") {
        # Minkowski distance (p=3 for demonstration)
        d <- dist(x, method = "minkowski", p = 3)
      } else {
        # Standard distance methods
        d <- dist(x, method = distance_method)
      }
      
      # Clustering
      if (requireNamespace("fastcluster", quietly = TRUE)) {
        hc <- fastcluster::hclust(d, method = "ward.D2")
      } else {
        hc <- hclust(d, method = "ward.D2")
      }
      
      return(hc)
      
    }, error = function(e) {
      warning(paste("Clustering with", distance_method, "failed:", e$message))
      # Fallback to euclidean
      d <- dist(x, method = "euclidean")
      return(hclust(d, method = "ward.D2"))
    })
  }
}

# Create clustering functions for different distance metrics
distance_methods <- c("euclidean", "manhattan", "canberra", "minkowski", "cosine")

print("Available distance metrics:")
for(method in distance_methods) {
  print(paste("-", method))
}

# ----------------------------
# Step 6: Create heatmaps with different distance metrics
# ----------------------------

# Choose which distance metric to use (change this line to test different metrics)
selected_distance <- "euclidean"  # Options: "euclidean", "manhattan", "canberra", "minkowski", "cosine"

print(paste("Creating heatmap with", selected_distance, "distance..."))

clustering_func <- create_clustering_function(selected_distance)

# Create individual heatmaps for each chromosome with PROMINENT BOXES
chr_heatmaps <- list()

for(chr_name in unique_chrs) {
  # Get segments for this chromosome
  chr_segments <- which(as.character(chromosomes) == chr_name)
  chr_matrix <- cna_matrix[chr_segments, , drop = FALSE]
  
  print(paste("Creating heatmap for chromosome", chr_name, "with", length(chr_segments), "segments"))
  
  # Create heatmap for this chromosome with STRONG BORDERS
  chr_heatmaps[[chr_name]] <- Heatmap(
    t(chr_matrix),
    name = if(chr_name == unique_chrs[1]) "Segment Mean" else NULL,
    col = col_fun,
    cluster_rows = if(chr_name == unique_chrs[1]) clustering_func else FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    # Enhanced chromosome labels
    column_title = chr_name,
    column_title_gp = gpar(fontsize = 16, fontface = "bold", col = "black"),
    column_title_side = "top",
    
    # Annotations
    top_annotation = NULL,
    left_annotation = if(chr_name == unique_chrs[1]) row_anno else NULL,
    
    # Visual settings
    use_raster = TRUE,
    raster_quality = 5,
    raster_device = "png",
    raster_device_param = list(type = "cairo", antialias = "default"),
    
    # Legend (only for first chromosome)
    show_heatmap_legend = if(chr_name == unique_chrs[1]) TRUE else FALSE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 12),
      border = "black",
      legend_height = unit(60, "mm"),
      legend_width = unit(8, "mm")
    ),
    
    # PROMINENT BORDERS - This creates the boxes around each chromosome
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 3),  # Thick black borders for clear boxes
    
    # Remove internal cell borders for cleaner look
    rect_gp = gpar(col = NA, lwd = 0)
  )
}

# Combine all chromosome heatmaps
print("Combining chromosome heatmaps...")
combined_ht <- chr_heatmaps[[1]]

if(length(chr_heatmaps) > 1) {
  for(i in 2:length(chr_heatmaps)) {
    combined_ht <- combined_ht + chr_heatmaps[[i]]
  }
}

# ----------------------------
# Step 7: Save the chromosome-boxed heatmap as PNG
# ----------------------------

output_file <- paste0("genome_wide_segm_means_heatmaps_per_chromosome_CLONE.", selected_distance, "_distance.png")

tryCatch({
  png(output_file, width = 6000, height = 2400, res = 300, type = "cairo", antialias = "default")
  
  print(paste("Drawing chromosome-separated heatmap with", selected_distance, "distance..."))
  
  # Add title
  grid.text("Genome-wide Copy Number Alterations by Chromosome (Clone Data)", 
            x = 0.5, y = 0.95, 
            gp = gpar(fontsize = 20, fontface = "bold"))
  
  grid.text(paste("Clustering Method:", toupper(selected_distance), "Distance"), 
            x = 0.5, y = 0.92, 
            gp = gpar(fontsize = 14, fontface = "italic", col = "gray30"))
  
  # Draw the heatmap with enhanced spacing between boxes
  draw(combined_ht, 
       heatmap_legend_side = "left", 
       annotation_legend_side = "left",
       gap = unit(4, "mm"),     # Larger gaps between chromosome boxes
       ht_gap = unit(4, "mm"),  # Additional spacing
       newpage = FALSE)
  
  dev.off()
  print(paste("Chromosome-boxed heatmap saved as:", output_file))
  
}, error = function(e) {
  if (dev.cur() != 1) dev.off()
  print(paste("Error creating chromosome-boxed heatmap:", e$message))
  
  # Fallback version with column splits
  print("Creating fallback version with chromosome boxes...")
  
  output_file_simple <- paste0("genome_wide_segm_means_heatmaps_per_chromosome_CLONE.", selected_distance, "_distance.fallback.png")
  png(output_file_simple, width = 5500, height = 2000, res = 300, type = "cairo", antialias = "default")
  
  ht_split <- Heatmap(
    t(cna_matrix),
    name = "Segment Mean",
    col = col_fun,
    cluster_rows = clustering_func,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    left_annotation = row_anno,
    use_raster = TRUE,
    raster_quality = 5,
    raster_device = "png",
    raster_device_param = list(type = "cairo", antialias = "default"),
    
    # Split by chromosome with STRONG BORDERS
    column_split = chromosomes,
    column_gap = unit(3, "mm"),  # Space between chromosome boxes
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 2),  # Thick borders around each chromosome
    
    column_title_side = "top",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 12),
      border = "black"
    ),
    
    # Clean appearance
    rect_gp = gpar(col = NA, lwd = 0)  # No internal cell borders
  )
  
  # Add title
  grid.text("Genome-wide Copy Number Alterations by Chromosome (Clone Data)", 
            x = 0.5, y = 0.95, 
            gp = gpar(fontsize = 18, fontface = "bold"))
  
  grid.text(paste("Clustering Method:", toupper(selected_distance), "Distance"), 
            x = 0.5, y = 0.92, 
            gp = gpar(fontsize = 12, fontface = "italic", col = "gray30"))
  
  draw(ht_split, 
       heatmap_legend_side = "left", 
       annotation_legend_side = "left",
       newpage = FALSE)
  
  dev.off()
  
  print(paste("Fallback chromosome-boxed heatmap saved as:", output_file_simple))
  output_file <- output_file_simple
})

# ----------------------------
# Optional: Create heatmaps for ALL distance metrics - ENABLED
# ----------------------------

create_all_distance_heatmaps <- function() {
  print("=== Creating heatmaps for all distance metrics ===")
  
  for(dist_method in distance_methods) {
    print(paste("Processing", dist_method, "distance..."))
    
    tryCatch({
      # Create clustering function for this distance
      clustering_func_temp <- create_clustering_function(dist_method)
      
      # Create filename
      output_file_temp <- paste0("genome_wide_segm_means_heatmaps_per_chromosome_CLONE.", dist_method, "_distance.png")
      
      # Create and save heatmap
      png(output_file_temp, width = 6000, height = 2400, res = 300, type = "cairo", antialias = "default")
      
      # Create individual heatmaps for each chromosome for this distance metric
      chr_heatmaps_temp <- list()
      
      for(chr_name in unique_chrs) {
        # Get segments for this chromosome
        chr_segments <- which(as.character(chromosomes) == chr_name)
        chr_matrix <- cna_matrix[chr_segments, , drop = FALSE]
        
        # Create heatmap for this chromosome
        chr_heatmaps_temp[[chr_name]] <- Heatmap(
          t(chr_matrix),
          name = if(chr_name == unique_chrs[1]) "Segment Mean" else NULL,
          col = col_fun,
          cluster_rows = if(chr_name == unique_chrs[1]) clustering_func_temp else FALSE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          
          # Enhanced chromosome labels
          column_title = chr_name,
          column_title_gp = gpar(fontsize = 16, fontface = "bold", col = "black"),
          column_title_side = "top",
          
          # Annotations
          top_annotation = NULL,
          left_annotation = if(chr_name == unique_chrs[1]) row_anno else NULL,
          
          # Visual settings
          use_raster = TRUE,
          raster_quality = 5,
          raster_device = "png",
          raster_device_param = list(type = "cairo", antialias = "default"),
          
          # Legend (only for first chromosome)
          show_heatmap_legend = if(chr_name == unique_chrs[1]) TRUE else FALSE,
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 14, fontface = "bold"),
            labels_gp = gpar(fontsize = 12),
            border = "black",
            legend_height = unit(60, "mm"),
            legend_width = unit(8, "mm")
          ),
          
          # PROMINENT BORDERS
          border = TRUE,
          border_gp = gpar(col = "black", lwd = 3),
          
          # Clean appearance
          rect_gp = gpar(col = NA, lwd = 0)
        )
      }
      
      # Combine all chromosome heatmaps for this distance
      combined_ht_temp <- chr_heatmaps_temp[[1]]
      if(length(chr_heatmaps_temp) > 1) {
        for(i in 2:length(chr_heatmaps_temp)) {
          combined_ht_temp <- combined_ht_temp + chr_heatmaps_temp[[i]]
        }
      }
      
      # Add title
      grid.text("Genome-wide Copy Number Alterations by Chromosome (Clone Data)", 
                x = 0.5, y = 0.95, 
                gp = gpar(fontsize = 20, fontface = "bold"))
      
      grid.text(paste("Clustering Method:", toupper(dist_method), "Distance"), 
                x = 0.5, y = 0.92, 
                gp = gpar(fontsize = 14, fontface = "italic", col = "gray30"))
      
      # Draw the heatmap
      draw(combined_ht_temp, 
           heatmap_legend_side = "left", 
           annotation_legend_side = "left",
           gap = unit(4, "mm"),
           ht_gap = unit(4, "mm"),
           newpage = FALSE)
      
      dev.off()
      
      print(paste("✓", dist_method, "distance heatmap saved as:", output_file_temp))
      
    }, error = function(e) {
      if (dev.cur() != 1) dev.off()
      print(paste("✗ Error with", dist_method, ":", e$message))
    })
  }
}

# AUTOMATICALLY CREATE ALL DISTANCE HEATMAPS
create_all_distance_heatmaps()

# ----------------------------
# Summary
# ----------------------------

print("=== Multi-Distance Chromosome-Boxed CNA Heatmap Complete (CLONE DATA) ===")

if (file.exists(output_file)) {
  file_size <- round(file.size(output_file) / 1024^2, 2)
  print(paste("✓ Heatmap created successfully:", output_file))
  print(paste("  Distance metric used:", selected_distance))
  print(paste("  File size:", file_size, "MB"))
  print(paste("  Format: High-quality PNG (300 DPI)"))
} else {
  stop("Error: Output file was not created")
}

print("=== Distance Metrics Available ===")
print("✓ Euclidean: Standard geometric distance")
print("✓ Manhattan: Sum of absolute differences (L1 norm)")
print("✓ Canberra: Weighted version of Manhattan distance")
print("✓ Minkowski: Generalized distance metric (p=3)")
print("✓ Cosine: Based on angle between vectors")

print("=== Usage Instructions ===")
print("To change distance metric, modify this line:")
print('selected_distance <- "euclidean"  # Change to: manhattan, canberra, minkowski, or cosine')
print("To create all distance versions, uncomment: create_all_distance_heatmaps()")

# print(paste("Data dimensions: cells =", ncol(cna_matrix), ", genomic segments =", nrow(cna_matrix)))
# print(paste("Chromosomes visualized:", paste(unique_chrs, collapse = ", ")))
# print(paste("Subclones:", paste(names(subclone_colors), collapse = ", ")))