# Step 1: Create DESeq2 Dataset
################################################################################
library(DESeq2)
library(EnhancedVolcano)
library(NMF)
library(dplyr)

# load in featureCounts file
readcounts_raw = read.table("featureCounts_P35.txt", header = TRUE, sep = "")


### clean readcounts into a numeric matrix with proper colnames
# set rownames to Gene ID and remove non-numeric columns
row.names(readcounts_raw) = readcounts_raw$Geneid
readcounts_numeric = readcounts_raw[ , -c(1:6)]

# extract SRR's from colnames, and rename cols using only SRR's
sample = gsub(".*STAR_2\\.([A-Za-z0-9]+)_Aligned.*", "\\1", names(readcounts_numeric))
colnames(readcounts_numeric) = sample

# create df of desired SRR's with corresponding cell types and replicate numbers
library(dplyr)
sample_df = tribble(
  ~SRR,             ~cell_type,       ~rep_num,
  "SRR10913364",    "EEC",            2,
  "SRR10913362",    "EEC",            1,
  "SRR10913360",    "Enterocyte",     5,
  "SRR10913358",    "Enterocyte",     4,
  "SRR10913356",    "Enterocyte",     3,
  "SRR10913354",    "Enterocyte",     2,
  "SRR10913352",    "Enterocyte",     1,
  "SRR10913349",    "Tuft",           5,
  "SRR10913346",    "Tuft",           4,
  "SRR10913343",    "Tuft",           3,
  "SRR10913340",    "Tuft",           2,
  "SRR10913337",    "Tuft",           1,
  "SRR10913334",    "SecPDG",         4,
  "SRR10913331",    "SecPDG",         3,
  "SRR10913328",    "SecPDG",         2,
  "SRR10913325",    "SecPDG",         1,
  "SRR10913323",    "AbsPro",         3,
  "SRR10913321",    "AbsPro",         2,
  "SRR10913319",    "AbsPro",         1,
  "SRR10913317",    "Stem",           3,
  "SRR10913315",    "Stem",           2,
  "SRR10913313",    "Stem",           1,
)

# drop unnecessary cols (duplicates with multiple runs)
readcounts_filtered = readcounts_numeric[, colnames(readcounts_numeric) %in% sample_df$SRR]

# reverse order of sample_df to match SRR's in readcounts table
flipped_sample_df <- sample_df[rev(rownames(sample_df)), ]

# conventionally name columns in readcounts table
readcounts = readcounts_filtered
names(readcounts) = paste(flipped_sample_df$cell_type, flipped_sample_df$rep_num, sep = ".")

# remove gene ID version numbers to avoid later problems
rownames(readcounts) <- sub("\\..*", "", rownames(readcounts))


### create the metadata
sample_info = data.frame(
  condition = flipped_sample_df$cell_type,
  row.names = names(readcounts))


### create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = readcounts,
  colData = sample_info,
  design = ~condition
)



# Step 2: Differential Expression Analysis
################################################################################
### Volcano Plot ###
# change reference level to compare Stem and Enterocyte
dds$condition = relevel(dds$condition, ref = "Stem")

# DESeq2 analysis
dds = DESeq(dds)

# Extract results for specific comparison
results_table <- results(dds, contrast = c("condition", "Stem", "Enterocyte"))


# Shrink log fold changes
resLFC <- lfcShrink(dds,
                    contrast = c("condition", "Stem", "Enterocyte"),
                    type = "normal")

# Generate basic volcano plot
EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = "log2FoldChange",
                y = "padj")

# Customized volcano plot with stricter thresholds
volcano_plot = EnhancedVolcano(resLFC,
                               lab = rownames(resLFC),
                               x = "log2FoldChange",
                               y = "padj",
                               pCutoff = 10e-45,
                               FCcutoff = 2)

# Save the plot to a file
ggsave("volcano_P35.pdf", plot = volcano_plot, width = 8, height = 6)

# pos log2FC: upregulation in Entrocytes


### Heatmap ###
# Perform rlog normalization
dds_rlog = rlog(dds, blind=FALSE)

# rlog-normalized counts
rlog_counts = assay(dds_rlog)

# Get differentially expressed genes
sigGenes = subset(resLFC, padj < 0.05)
DEgenes = rownames(sigGenes)

# Get normalized counts for differentially expressed genes
mat_genes_all = rlog_counts[DEgenes, c('Stem.1', 'Stem.2', 'Stem.3', 
                                       'Enterocyte.1', 'Enterocyte.2', 
                                       'Enterocyte.3', 'Enterocyte.4', 
                                       'Enterocyte.5')]

# Select top 200 most variable DE genes
vars = apply(mat_genes_all, 1, var)
top200_idx = order(vars, decreasing = TRUE)[1:200]
mat_genes = mat_genes_all[top200_idx, ]

# Create the heatmap with correlation-based clustering
aheatmap(mat_genes, Rowv = TRUE, Colv = NA,
         distfun = "correlation", hclustfun = "complete",
         width = 8.5, height = 11, scale = "row",
         filename = "heatmap_P35.pdf")


# Step 3: create files for functional analysis
################################################################################
### Generate the .gct file
# Extract gene ID's to convert to gene symbols
gene_list <- rownames(rlog_counts)
writeLines(gene_list, "genes_for_gprofiler_P35.txt")

# Import the geneID-symbol conversion table
conversion_table <- read.csv("converted_genes_P35.csv")

# Merge the rlog-transformed data with the geneID-symbol table
merged_data <- merge(conversion_table, rlog_counts, 
                     by.x = "initial_alias", by.y = "row.names")

# Replace the geneID column with gene symbols, and remove all other non-data columns
final_data <- merged_data[, colnames(merged_data)[-c(1, 2, 4, 5)]]

# Convert the gene symbols to uppercase
# final_data$name <- toupper(final_data$name)

# Add a 'description' column (can be blank or "NA")
final_data <- final_data %>%
  mutate(description = "NA") %>%
  relocate(name, description)

# Check dimensions
cat("Rows:", nrow(final_data), "Columns (excluding name/desc):", ncol(final_data) - 2, "\n")

# ---- Custom function to write a valid GCT file ----
write_gct <- function(df, file) {
  con <- file(file, "w")
  
  # Write GCT header lines
  writeLines("#1.2", con)
  writeLines(sprintf("%d\t%d", nrow(df), ncol(df) - 2), con)
  
  # Write data table
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  close(con)
}

# Write the final .gct file
write_gct(final_data, "output_P35.gct")


### Generate the .cls file
# First line: <num_samples> <num_classes> 1
row1_cls <- "22 6 1"

# Second line: class names, preceded by "#"
row2_cls <- "# Stem AbsPro SecPDG Tuft Enterocyte EEC"

# Third line: group labels (indices for each sample)
row3_cls <- "0 0 0 1 1 1 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5"

# Combine all rows into one character vector
cls_final <- c(row1_cls, row2_cls, row3_cls)

# Write to file line-by-line
writeLines(cls_final, "output_P35.cls")