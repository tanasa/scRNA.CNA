# Fixed Clean CNA Heatmap Script with Multiple Distance Metrics
# revised by Claude Sonnet 

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# Check if required objects exist
required_objects <- c("CNA_mtx_relat_annot", "obj_scevan2")
missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
  stop(paste("Missing required objects:", paste(missing_objects, collapse = ", ")))
}

# ----------------------------
# Step 1: Prepare CNA matrix
# ----------------------------

# Check if CNA_mtx_relat_annot exists and has the right structure
if (!exists("CNA_mtx_relat_annot")) {
  # If it doesn't exist, create it from the existing objects
  if (exists("CNA_mtx_relat") && exists("count_mtx_annot")) {
    cat("Creating CNA_mtx_relat_annot from existing objects...\n")
    CNA_mtx_relat_annot <- cbind(count_mtx_annot, CNA_mtx_relat)
  } else {
    stop("Required objects CNA_mtx_relat_annot, CNA_mtx_relat, or count_mtx_annot not found")
  }
}

# Extract genomic annotation columns (first 5 columns)
genomic_info <- CNA_mtx_relat_annot[, 1:5]

# Extract CNA values (columns 6 onwards - the cell data)
cna_matrix <- CNA_mtx_relat_annot[, 6:ncol(CNA_mtx_relat_annot)]

# Get metadata - fix the metadata reference
if (exists("obj_scevan2") && "meta.data" %in% slotNames(obj_scevan2)) {
  metadata <- obj_scevan2@meta.data
} else if (exists("metadata")) {
  # metadata already exists
  metadata <- metadata
} else {
  stop("Cannot find metadata. Expected obj_scevan2@meta.data or metadata object")
}

# Ensure we have matching cells between CNA matrix and metadata
common_cells <- intersect(colnames(cna_matrix), rownames(metadata))

if (length(common_cells) == 0) {
  stop("No common cells found between CNA matrix and metadata")
}

cna_matrix <- cna_matrix[, common_cells, drop = FALSE]
metadata_filtered <- metadata[common_cells, , drop = FALSE]

cat("Number of cells:", ncol(cna_matrix), "\n")
cat("Number of genes/segments:", nrow(cna_matrix), "\n")

# ----------------------------
# Step 2: Prepare chromosome separation and colors
# ----------------------------

# Ensure genomic_info has the right column names
if (!"seqnames" %in% colnames(genomic_info)) {
  if ("chr" %in% colnames(genomic_info)) {
    colnames(genomic_info)[colnames(genomic_info) == "chr"] <- "seqnames"
  } else if (ncol(genomic_info) >= 1) {
    colnames(genomic_info)[1] <- "seqnames"
  } else {
    stop("Cannot identify chromosome column in genomic_info")
  }
}

# Create chromosome factor with proper ordering
chr_levels <- c(as.character(1:22), "X", "Y")
chromosomes <- factor(genomic_info$seqnames, levels = chr_levels)

# Check what chromosomes are actually present in your data
unique_chrs <- unique(as.character(chromosomes))
unique_chrs <- unique_chrs[!is.na(unique_chrs)]  # Remove NA values

cat("Chromosomes in data:", paste(unique_chrs, collapse = ", "), "\n")

if (length(unique_chrs) == 0) {
  stop("No valid chromosomes found in genomic_info")
}

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

cat("Chromosome boundaries calculated\n")

# ----------------------------
# Step 3: Prepare row annotation (cell classification)
# ----------------------------

# Extract cell classifications - handle missing class column
if ("class" %in% colnames(metadata_filtered)) {
  cell_classes <- metadata_filtered$class
} else {
  cat("Warning: 'class' column not found in metadata. Using dummy classification.\n")
  cell_classes <- rep("unknown", nrow(metadata_filtered))
}

names(cell_classes) <- rownames(metadata_filtered)

# Check what classes are present
unique_classes <- unique(cell_classes[!is.na(cell_classes)])
cat("Cell classes in data:", paste(unique_classes, collapse = ", "), "\n")

# Create color palette for cell classes
if(length(unique_classes) <= 8 && length(unique_classes) > 0) {
  class_colors <- RColorBrewer::brewer.pal(max(3, length(unique_classes)), "Dark2")[1:length(unique_classes)]
} else if (length(unique_classes) > 8) {
  class_colors <- rainbow(length(unique_classes))
} else {
  class_colors <- "gray"
}
names(class_colors) <- unique_classes

# Handle NA values
if(any(is.na(cell_classes))) {
  class_colors <- c(class_colors, "white")
  names(class_colors)[length(class_colors)] <- "NA"
  cell_classes[is.na(cell_classes)] <- "NA"
}

cat("Class colors assigned\n")

# Create row annotation
tryCatch({
  row_anno <- rowAnnotation(
    Class = cell_classes,
    col = list(Class = class_colors),
    annotation_name_side = "top",
    border = TRUE,
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
  )
}, error = function(e) {
  cat("Error creating row annotation:", e$message, "\n")
  row_anno <- NULL
})

# ----------------------------
# Step 4: Prepare color scale for CNA values
# ----------------------------

# Calculate data range - convert to numeric matrix safely
cna_values <- as.matrix(cna_matrix)
mode(cna_values) <- "numeric"  # Ensure numeric

# Remove any infinite or non-finite values
cna_values[!is.finite(cna_values)] <- NA

data_range <- range(cna_values, na.rm = TRUE)
cat("Data range:", round(data_range[1], 3), "to", round(data_range[2], 3), "\n")

# Fixed color scale between -1 to +1
col_breaks <- c(-1, -0.5, -0.2, 0, 0.2, 0.5, 1)

tryCatch({
  col_fun <- colorRamp2(col_breaks, 
                        colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(n = length(col_breaks)))
}, error = function(e) {
  cat("Error creating color function, using simple gradient\n")
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
})

cat("Color scale: -1 to +1\n")

# ----------------------------
# Step 5: Define multiple distance metrics clustering functions
# ----------------------------

cat("Setting up multiple distance metrics...\n")

# Function to compute cosine distance safely
cosine_distance <- function(x, y) {
  tryCatch({
    dot_product <- sum(x * y, na.rm = TRUE)
    norm_x <- sqrt(sum(x^2, na.rm = TRUE))
    norm_y <- sqrt(sum(y^2, na.rm = TRUE))
    
    if (norm_x == 0 || norm_y == 0) {
      return(1)  # Maximum distance for zero vectors
    }
    
    cosine_sim <- dot_product / (norm_x * norm_y)
    return(1 - cosine_sim)
  }, error = function(e) {
    return(1)  # Return maximum distance on error
  })
}

# Advanced clustering function with multiple distance metrics
create_clustering_function <- function(distance_method = "euclidean") {
  function(x) {
    tryCatch({
      # Remove rows/columns with all NA
      complete_rows <- complete.cases(x)
      if (sum(complete_rows) < 2) {
        stop("Too few complete cases for clustering")
      }
      
      x_clean <- x[complete_rows, , drop = FALSE]
      
      # Different distance calculations
      if (distance_method == "cosine") {
        # Custom cosine distance implementation
        n <- nrow(x_clean)
        d_matrix <- matrix(0, n, n)
        for(i in 1:(n-1)) {
          for(j in (i+1):n) {
            d_matrix[i,j] <- d_matrix[j,i] <- cosine_distance(x_clean[i,], x_clean[j,])
          }
        }
        d <- as.dist(d_matrix)
      } else if (distance_method == "minkowski") {
        # Minkowski distance (p=3 for demonstration)
        d <- dist(x_clean, method = "minkowski", p = 3)
      } else {
        # Standard distance methods
        d <- dist(x_clean, method = distance_method)
      }
      
      # Clustering
      if (requireNamespace("fastcluster", quietly = TRUE)) {
        hc <- fastcluster::hclust(d, method = "ward.D2")
      } else {
        hc <- hclust(d, method = "ward.D2")
      }
      
      # Create full result with correct order
      full_hc <- hc
      if (sum(complete_rows) < nrow(x)) {
        # Handle incomplete cases by creating a full dendrogram
        full_order <- rep(NA, nrow(x))
        full_order[complete_rows] <- hc$order
        incomplete_indices <- which(!complete_rows)
        full_order[incomplete_indices] <- (max(hc$order, na.rm = TRUE) + 1):(max(hc$order, na.rm = TRUE) + length(incomplete_indices))
        full_hc$order <- full_order[!is.na(full_order)]
      }
      
      return(full_hc)
      
    }, error = function(e) {
      warning(paste("Clustering with", distance_method, "failed:", e$message))
      # Fallback to simple ordering
      n <- nrow(x)
      fake_hc <- list(
        merge = matrix(c(-1, -2), nrow = 1),
        height = 1,
        order = 1:n,
        labels = rownames(x),
        call = call("hclust"),
        method = "ward.D2",
        dist.method = "euclidean"
      )
      class(fake_hc) <- "hclust"
      return(fake_hc)
    })
  }
}

# Create clustering functions for different distance metrics
distance_methods <- c("euclidean", "manhattan", "canberra", "minkowski", "cosine")

cat("Available distance metrics:\n")
for(method in distance_methods) {
  cat("-", method, "\n")
}

# ----------------------------
# Step 6: Create heatmaps with different distance metrics
# ----------------------------

# Choose which distance metric to use (change this line to test different metrics)
selected_distance <- "euclidean"  # Options: "euclidean", "manhattan", "canberra", "minkowski", "cosine"

cat("Creating heatmap with", selected_distance, "distance...\n")

clustering_func <- create_clustering_function(selected_distance)

# Create individual heatmaps for each chromosome with PROMINENT BOXES
chr_heatmaps <- list()

for(chr_name in unique_chrs) {
  # Get segments for this chromosome
  chr_segments <- which(as.character(chromosomes) == chr_name)
  chr_matrix <- cna_matrix[chr_segments, , drop = FALSE]
  
  cat("Creating heatmap for chromosome", chr_name, "with", length(chr_segments), "segments\n")
  
  # Create heatmap for this chromosome with STRONG BORDERS
  tryCatch({
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
      left_annotation = if(chr_name == unique_chrs[1] && !is.null(row_anno)) row_anno else NULL,
      
      # Visual settings
      use_raster = TRUE,
      raster_quality = 5,
      
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
  }, error = function(e) {
    cat("Error creating heatmap for chromosome", chr_name, ":", e$message, "\n")
    # Create a simple placeholder
    chr_heatmaps[[chr_name]] <- NULL
  })
}

# Remove NULL entries
chr_heatmaps <- chr_heatmaps[!sapply(chr_heatmaps, is.null)]

if (length(chr_heatmaps) == 0) {
  stop("No heatmaps could be created")
}

# Combine all chromosome heatmaps
cat("Combining chromosome heatmaps...\n")
combined_ht <- chr_heatmaps[[1]]

if(length(chr_heatmaps) > 1) {
  for(i in 2:length(chr_heatmaps)) {
    combined_ht <- combined_ht + chr_heatmaps[[i]]
  }
}

# ----------------------------
# Step 7: Save the chromosome-boxed heatmap as PNG
# ----------------------------

output_file <- paste0("genome_wide_segm_means_heatmaps_per_chromosome.", selected_distance, "_distance.png")

tryCatch({
  png(output_file, width = 6000, height = 2400, res = 300, type = "cairo", antialias = "default")
  
  cat("Drawing chromosome-separated heatmap with", selected_distance, "distance...\n")
  
  # Add title
  grid.text("Genome-wide Copy Number Alterations by Chromosome", 
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
  cat("Chromosome-boxed heatmap saved as:", output_file, "\n")
  
}, error = function(e) {
  if (dev.cur() != 1) dev.off()
  cat("Error creating chromosome-boxed heatmap:", e$message, "\n")
  
  # Fallback version with column splits
  cat("Creating fallback version with chromosome boxes...\n")
  
  output_file_simple <- paste0("genome_wide_segm_means_heatmaps_per_chromosome.", selected_distance, "_distance.fallback.png")
  
  tryCatch({
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
    grid.text("Genome-wide Copy Number Alterations by Chromosome", 
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
    
    cat("Fallback chromosome-boxed heatmap saved as:", output_file_simple, "\n")
    output_file <- output_file_simple
    
  }, error = function(e2) {
    if (dev.cur() != 1) dev.off()
    cat("Both heatmap creation methods failed:", e2$message, "\n")
  })
})

# ----------------------------
# Summary
# ----------------------------

cat("=== Multi-Distance Chromosome-Boxed CNA Heatmap Complete ===\n")

if (exists("output_file") && file.exists(output_file)) {
  file_size <- round(file.size(output_file) / 1024^2, 2)
  cat("✓ Heatmap created successfully:", output_file, "\n")
  cat("  Distance metric used:", selected_distance, "\n")
  cat("  File size:", file_size, "MB\n")
  cat("  Format: High-quality PNG (300 DPI)\n")
} else {
  cat("✗ Error: Output file was not created\n")
}

cat("Data dimensions: cells =", ncol(cna_matrix), ", genomic segments =", nrow(cna_matrix), "\n")
cat("Chromosomes visualized:", paste(unique_chrs, collapse = ", "), "\n")
if (exists("unique_classes")) {
  cat("Cell classes:", paste(unique_classes, collapse = ", "), "\n")
}
