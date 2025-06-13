# Load libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# --- Step 1: Transpose CNA matrix ---
mat <- t(as.matrix(oncoHeat))  # oncoHeat: rows = clones, cols = bands

# --- Step 2: Automatically define subclone_vector ---
subclone_vector <- paste0("Clone", seq_len(nrow(oncoHeat)))  # e.g., "Clone1", "Clone2", ...
stopifnot(length(subclone_vector) == ncol(mat))

# --- Step 3: Generate dynamic clone colors ---
subclone_levels <- unique(subclone_vector)
n_clones <- length(subclone_levels)

if (n_clones <= 12) {
  palette_base <- brewer.pal(n = max(n_clones, 3), name = "Set3")
} else {
  palette_base <- colorRampPalette(brewer.pal(12, "Set3"))(n_clones)
}
subclone_colors <- setNames(palette_base[1:n_clones], subclone_levels)

# --- Step 4: Define CNA color scale ---
col_fun <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("blue", "lightblue", "white", "pink", "red")
)

# --- Step 5: Subclone annotation bar ---
ha <- HeatmapAnnotation(
  Subclone = subclone_vector,
  col = list(Subclone = subclone_colors),
  annotation_name_side = "left",
  show_annotation_name = TRUE
)

# --- Step 6: Plot heatmap ---
ht <- Heatmap(
  mat,
  name = "CNA",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = ha,
  row_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title = "Copy Number",
    at = c(-2, -1, 0, 1, 2),
    labels = c("DEL", "LOSS", "Neutral", "GAIN", "AMP")
  )
)

# --- Step 7: Export output with taller height and narrower width ---

# PNG: width reduced from 2000 → 1200 pixels, height increased to 2000
png("oncoHeatmap_cytobands.png", width = 800, height = 2000, res = 200)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# PDF: width reduced from 14 → 8 inches, height increased to 12
pdf("oncoHeatmap_cytobands.pdf", width = 6, height = 12)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

