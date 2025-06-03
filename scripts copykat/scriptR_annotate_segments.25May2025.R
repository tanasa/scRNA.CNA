#!/usr/bin/env Rscript

# === Load libraries ===
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(AnnotationDbi)
library(AnnotationHub)
library(GenomeInfoDb)

# === 1. Load CNA segments ===
# args <- commandArgs(trailingOnly = TRUE)

# if (length(args) == 0) {
#    stop("Please provide the CNA segments file as an argument.")
# }

# seg_file <- args[1]
seg_file="LGG_copykat_CNA_results.seg"

cat("Reading CNA segments from:", seg_file, "\n")

seg <- read.table(seg_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Standardize chromosome naming to UCSC style
seg$chrom <- ifelse(grepl("^chr", seg$chrom), seg$chrom, paste0("chr", seg$chrom))
# Fix common alias for chrX
seg$chrom[seg$chrom == "chr23"] <- "chrX"

# Filter by log2 ratio threshold
seg_filtered <- seg[seg$seg.mean > 0.2 | seg$seg.mean < -0.2, ]

# Save filtered segments
write.table(seg_filtered, 
            "copykat_CNA_filtered_segments_log2R_0.2.tsv",
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Convert to GRanges
seg_gr <- GRanges(
  seqnames = seg_filtered$chrom,
  ranges = IRanges(start = seg_filtered$loc.start, end = seg_filtered$loc.end),
  ID = seg_filtered$ID,
  seg.mean = seg_filtered$seg.mean
)

# === 2. Get gene coordinates from TxDb ===
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)

# Add gene symbols using org.Hs.eg.db

gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = as.character(genes$gene_id),
                                      columns = "SYMBOL",
                                      keytype = "ENTREZID")

mcols(genes)$SYMBOL <- gene_symbols$SYMBOL[match(as.character(genes$gene_id),
                                                  gene_symbols$ENTREZID)]

# === 3. Harmonize seqlevels between segments and genes ===
seqlevelsStyle(seg_gr) <- "UCSC"
seqlevelsStyle(genes) <- "UCSC"

# Remove unused sequence levels
seqlevels(seg_gr, pruning.mode = "coarse") <- seqlevelsInUse(seg_gr)
seqlevels(genes, pruning.mode = "coarse") <- seqlevelsInUse(genes)

# === 4. Find overlaps between CNA segments and genes ===
ov <- findOverlaps(seg_gr, genes)
cat("Number of overlaps found: ", length(ov), "\n")

# === 5. Conditional creation of the annotated segments dataframe ===

  annotated_segments <- data.frame(
    segment_ID  = seg_gr[queryHits(ov)]$ID,
    chrom       = as.character(seqnames(seg_gr[queryHits(ov)])),
    start       = start(seg_gr[queryHits(ov)]),
    end         = end(seg_gr[queryHits(ov)]),
    seg_mean    = seg_gr[queryHits(ov)]$seg.mean,
    gene_id     = genes[subjectHits(ov)]$gene_id,
    gene_symbol = genes[subjectHits(ov)]$SYMBOL,
    gene_start  = start(genes[subjectHits(ov)]),
    gene_end    = end(genes[subjectHits(ov)])
  )

  # === 6. Annotate with cytoband information ===
  ah <- AnnotationHub()
  cyto <- query(ah, c("hg38", "cytoBand"))[[1]]
  
  # Create GRanges object for genes
  gene_gr <- GRanges(
    seqnames = annotated_segments$chrom,
    ranges = IRanges(start = annotated_segments$gene_start,
                      end = annotated_segments$gene_end)
  )
  
  cyto_ov <- findOverlaps(gene_gr, cyto)
  annotated_segments$cytoband <- NA
  annotated_segments$cytoband[queryHits(cyto_ov)] <- mcols(cyto[subjectHits(cyto_ov)])$name
  
  # Annotate gain/loss type
  annotated_segments$type <- ifelse(annotated_segments$seg_mean > 0.2, "gain", "loss")
  
  head(annotated_segments, 3)
#                 segment_ID chrom    start      end   seg_mean   gene_id
# 1 LGG.03_AAACCCAAGCCTTTGA.1  chr1 37873972 42613844 -0.2593671 100130557
# 2 LGG.03_AAACCCAAGCCTTTGA.1  chr1 37873972 42613844 -0.2593671 100500801
# 3 LGG.03_AAACCCAAGCCTTTGA.1  chr1 37873972 42613844 -0.2593671 100507178
#  gene_symbol gene_start gene_end cytoband type
# 1    NFYC-AS1   40690380 40692074       85 loss
# 2     MIR3659   38089231 38089329       86 loss
# 3  SLFNL1-AS1   41014590 41043890       85 loss
  dim(annotated_segments)

  # Save full annotated table
  write.table(annotated_segments,
               "copykat_CNA_filtered_segments_log2R_0.2_annotated_segments_cytoband.tsv",
               sep = "\t", 
               row.names = FALSE, 
               quote = FALSE)
 
# a data frame where each row represents a gene overlapping the same segment 
# (for example, segment LGG.03_AAACCCAAGCCTTTGA.1 on chr1).
# You’d like to collapse this into a structure where each segment 
# has a list of overlapping genes.
 
library(dplyr)

# Collapse gene symbols per segment
collapsed_segments <- annotated_segments %>%
  group_by(segment_ID, chrom, start, end, seg_mean, type) %>%
  summarise(
    genes = paste(unique(gene_symbol), collapse = ", "),
    cytobands = paste(unique(cytoband), collapse = ", "),
    n_genes = n()
  ) %>%
  ungroup()

dim(collapsed_segments)

# head(collapsed_segments)
# A tibble: 6 × 8
#  segment_ID                chrom    start      end seg_mean type  genes n_genes
#  <chr>                     <chr>    <int>    <int>    <dbl> <chr> <chr>   <int>
# 1 LGG.03_AAACCCAAGCCTTTGA.1 chr1  37873972 42613844   -0.259 loss  NFYC…      76
# 2 LGG.03_AAACCCAAGCCTTTGA.1 chr1  42841006 43066488   -0.201 loss  ERMA…       5

  # Save full collapsed table
  write.table(collapsed_segments,
               "copykat_CNA_filtered_segments_log2R_0.2_collapsed_segments_cytoband.tsv",
               sep = "\t", 
               row.names = FALSE, 
               quote = FALSE)

  # === 7. Filter for known genes of interest (oncogenes & TSGs) ===
  
oncogenes <- c(
  "IDH1", "IDH2", "EGFR", "PIK3CA", "PIK3R1", "PDGFRA", "MYC", "MYCN",
  "MET", "CDK4", "MDM2", "MDM4", "TERT", "SOX2", "BRAF", "H3F3A",
  "NOTCH1", "NOTCH2", "GNAS", 
  # Additional oncogenes frequently seen in gliomas/GBM
  "FGFR1", "FGFR3", "KRAS", "NRAS", "HRAS", "AKT3", "CCND1"
)


tumor_suppressors <- c(
  "TP53", "ATRX", "CDKN2A", "CDKN2B", "PTEN", "NF1", "RB1", "TET2",
  "FUBP1", "CIC", "NDRG2", "GATA4", "ZEB1",
  # Additional TSGs frequently seen in gliomas/GBM
  "PBRM1", "SETD2", "ARID1A", "ARID1B", "SMARCB1", "SMARCA4", "LZTR1"
)

  genes_of_interest <- c(oncogenes, tumor_suppressors)

  matches <- annotated_segments %>%
    filter(gene_symbol %in% genes_of_interest)
  
  # Keep the most extreme log2 ratio per gene
  matches_gene_interest <- matches %>%
    group_by(gene_symbol) %>%
    slice_max(abs(seg_mean), with_ties = FALSE) %>%
    ungroup()
  
  # Save to file
  write.table(matches_gene_interest,
              "copykat_CNA_filtered_segments_log2R_0.2_genes_of_interest.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)


# head(annotated_segments,3)
#                 segment_ID chrom    start      end   seg_mean   gene_id
# 1 LGG.03_AAACCCAAGCCTTTGA.1  chr1 37873972 42613844 -0.2593671 100130557
# 2 LGG.03_AAACCCAAGCCTTTGA.1  chr1 37873972 42613844 -0.2593671 100500801
# 3 LGG.03_AAACCCAAGCCTTTGA.1  chr1 37873972 42613844 -0.2593671 100507178
#  gene_symbol gene_start gene_end cytoband type
# 1    NFYC-AS1   40690380 40692074       85 loss
# 2     MIR3659   38089231 38089329       86 loss
# 3  SLFNL1-AS1   41014590 41043890       85 loss
# head(collapsed_segments,3)
# A tibble: 3 × 8
#  segment_ID                chrom    start      end seg_mean type  genes n_genes
#  <chr>                     <chr>    <int>    <int>    <dbl> <chr> <chr>   <int>
# 1 LGG.03_AAACCCAAGCCTTTGA.1 chr1  37873972 42613844   -0.259 loss  NFYC…      76
# 2 LGG.03_AAACCCAAGCCTTTGA.1 chr1  42841006 43066488   -0.201 loss  ERMA…       5
# 3 LGG.03_AAACCCAAGCCTTTGA.1 chr12  7107665 11169762   -0.262 loss  MIR1…      92
