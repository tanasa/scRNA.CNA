# Suppress ComplexHeatmap messages
library(ComplexHeatmap) 
ht_opt$message = FALSE

# Blue to red color scale (better for subtle differences)
p_blue_red <- FeaturePlot(
  obj_subset,
  features = paste0(gene_of_interest, "_CNA"),
  pt.size = 0.5  # Reduced from 0.8 to 0.5
) + 
  scale_colour_distiller(
    palette = "RdBu", 
    direction = -1,  # -1 makes blue = low (negative), red = high (positive)
    name = "Segm Mean"
  ) +
  ggtitle(paste("Segm Mean")) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

print(p_blue_red)

# Save the plot with smaller dimensions
ggsave(
  filename = paste0("gene_", gene_of_interest, "_segm_mean_blue_red.png"),
  plot = p_blue_red,
  width = 6,   # Reduced from 10 to 6
  height = 5,  # Reduced from 8 to 5
  dpi = 300,
  bg = "white"
)

cat("Plot saved as:", paste0("gene_", gene_of_interest, "_segm_mean_blue_red.png"), "\n")

######################################################################################
######################################################################################

p_rdylbu <- FeaturePlot(
  obj_subset,
  features = paste0(gene_of_interest, "_CNA"),
  pt.size = 0.5  # Reduced from 0.8 to 0.5
) + 
  scale_colour_distiller(
    palette = "RdYlBu", 
    direction = -1,  # Blue = deletions, Red = amplifications
    name = "Segm Mean"
  ) +
  ggtitle(paste("Segm Mean")) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

print(p_rdylbu)

# Save the plot with smaller dimensions
ggsave(
  filename = paste0("gene_", gene_of_interest, "_segm_mean_red_yellow_blue.png"),
  plot = p_rdylbu,
  width = 6,   # Reduced from 10 to 6
  height = 5,  # Reduced from 8 to 5
  dpi = 300,
  bg = "white"
)

cat("Plot saved as:", paste0("gene_", gene_of_interest, "_segm_mean_red_yellow_blue.png"), "\n")

######################################################################################
######################################################################################

p_spectral <- FeaturePlot(
  obj_subset,
  features = paste0(gene_of_interest, "_CNA"),
  pt.size = 0.5  # Reduced from 0.8 to 0.5
) + 
  scale_colour_distiller(
    palette = "Spectral", 
    direction = -1,  # Blue = deletions, Red = amplifications
    name = "Segm Mean"
  ) +
  ggtitle(paste("Segm Mean")) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

print(p_spectral)

# Save the plot with smaller dimensions
ggsave(
  filename = paste0("gene_", gene_of_interest, "_segm_mean_spectral.png"),
  plot = p_spectral,
  width = 6,   # Reduced from 10 to 6
  height = 5,  # Reduced from 8 to 5
  dpi = 300,
  bg = "white"
)

cat("Plot saved as:", paste0("gene_", gene_of_interest, "_segm_mean_spectral.png"), "\n")

######################################################################################
######################################################################################

p_custom <- FeaturePlot(
  obj_subset,
  features = paste0(gene_of_interest, "_CNA"),
  pt.size = 0.5  # Reduced from 0.8 to 0.5
) + 
  scale_colour_gradientn(
    colors = c("darkblue", "blue", "lightblue", "white", "lightcoral", "red", "darkred"),
    name = "Segm Mean"
  ) +
  ggtitle(paste("Segm Mean")) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

print(p_custom)

# Save the plot with smaller dimensions
ggsave(
  filename = paste0("gene_", gene_of_interest, "_segm_mean_custom.png"),
  plot = p_custom,
  width = 6,   # Reduced from 10 to 6
  height = 5,  # Reduced from 8 to 5
  dpi = 300,
  bg = "white"
)

cat("Plot saved as:", paste0("gene_", gene_of_interest, "_segm_mean_custom.png"), "\n")

######################################################################################
######################################################################################

p_custom2 <- FeaturePlot(
  obj_subset,
  features = paste0(gene_of_interest, "_CNA"),
  pt.size = 0.5  # Reduced from 0.8 to 0.5
) + 
  scale_colour_gradientn(
    colors = c("darkblue", "blue", "lightblue", "white", "lightcoral", "red", "darkred"),
    values = scales::rescale(c(-1, -0.5, -0.2, 0, 0.2, 0.5, 1)),  # Updated rescale values
    limits = c(-1, 1),  # Set explicit limits from -1 to +1
    name = "Segm Mean"
  ) +
  ggtitle(paste("Segm Mean")) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

print(p_custom2)

# Save the plot with smaller dimensions
ggsave(
  filename = paste0("gene_", gene_of_interest, "_segm_mean_custom_scaled.png"),
  plot = p_custom2,
  width = 6,   # Reduced from 10 to 6
  height = 5,  # Reduced from 8 to 5
  dpi = 300,
  bg = "white"
)

cat("Plot saved as:", paste0("gene_", gene_of_interest, "_segm_mean_custom_scaled.png"), "\n")
