
# to generate gene positions file
library(rtracklayer)
library(data.table)
library(stringr)
library(dplyr)

# Step 1: Read the GTF file
gtf_file <- "hg38.refGene.gtf"
gtf <- fread(gtf_file, sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE)

# Step 2: Assign GTF column names
colnames(gtf) <- c("chr", "source", "feature", "start", "end", 
                   "score", "strand", "frame", "attributes")

# Step 3: Filter for transcript entries (gene-level proxies)
transcripts <- gtf %>% filter(feature == "transcript")

# Step 4: Extract gene_name from attributes
transcripts$gene_name <- str_match(transcripts$attributes, 'gene_name "([^"]+)"')[,2]

# Step 5: Build gene-level data frame
gene_df <- transcripts %>%
  select(Gene = gene_name, Chromosome = chr, Start = start, End = end)

# Step 6: Remove duplicate coordinate entries (same chr, start, end)
gene_df_unique <- gene_df %>%
  distinct(Chromosome, Start, End, .keep_all = TRUE)

# Step 7: Keep the longest isoform per gene
gene_df_unique$Length <- gene_df_unique$End - gene_df_unique$Start
gene_longest <- gene_df_unique %>%
  group_by(Gene) %>%
  slice_max(order_by = Length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-Length) %>%
  arrange(Chromosome, Start)

# Step 8: Write to output file
write.table(gene_longest, file = "hg38_gene_positions_longest_isoform.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("✔ Final gene position file (longest isoform per gene) saved to 'gene_positions_longest_isoform.tsv'\n")