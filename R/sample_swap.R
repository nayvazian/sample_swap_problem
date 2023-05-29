library(DESeq2)

readSample <- function(path) {
  counts <- read.table(path, header = FALSE, row.names = 1, sep = "\t", check.names = FALSE)
  counts <- tail(counts, -4)[,1, drop=FALSE]
  return(counts)
}

counts_df <- data.frame(readSample("C:/Users/Nick/gitbash_documents/sample_swap_problem/gene_expression_data/donor_A/colon.tsv"),
                        readSample("C:/Users/Nick/gitbash_documents/sample_swap_problem/gene_expression_data/donor_A/lung.tsv")$V2,
                        readSample("C:/Users/Nick/gitbash_documents/sample_swap_problem/gene_expression_data/donor_A/spleen.tsv")$V2,
                        readSample("C:/Users/Nick/gitbash_documents/sample_swap_problem/gene_expression_data/donor_A/stomach.tsv")$V2,
                        readSample("C:/Users/Nick/gitbash_documents/sample_swap_problem/gene_expression_data/donor_B/colon.tsv")$V2,
                        readSample("C:/Users/Nick/gitbash_documents/sample_swap_problem/gene_expression_data/donor_B/lung.tsv")$V2,
                        readSample("C:/Users/Nick/gitbash_documents/sample_swap_problem/gene_expression_data/donor_B/spleen.tsv")$V2,
                        readSample("C:/Users/Nick/gitbash_documents/sample_swap_problem/gene_expression_data/donor_B/stomach.tsv")$V2)
colnames(counts_df) <- c("donorA_colon", "donorA_lung", "donorA_spleen", "donorA_stomach", "donorB_colon", "donorB_lung", "donorB_spleen", "donorB_stomach")

metadata <- data.frame(samples=colnames(counts_df), 
                       donor=as.factor(c(rep('Donor_A',4), rep('Donor_B', 4))),
                       tissue=as.factor(c("colon","lung","spleen","stomach","colon","lung","spleen","stomach")), row.names = TRUE)

dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = metadata, design = ~ tissue + donor)
dds <- DESeq(dds)

results <- results(dds)
log2fc <- log2(results$baseMean + 1)  # Adding 1 to avoid log(0)
pvalues <- results$pvalue
padj <- results$padj

# Create a dataframe
results_df <- data.frame(log2fc, pvalues, padj)


Donor1_genes <- subset(results, padj < 0.05 & log2FoldChange > 0 & colData(dds)$Donor == "Donor1")
Donor2_genes <- subset(results, padj < 0.05 & log2FoldChange < 0 & colData(dds)$Donor == "Donor2")

# Find genes showing opposite LFCs between patients
swapped_genes <- intersect(row.names(Donor1_genes), row.names(Donor2_genes))

normalized_counts <- counts(dds, normalized = TRUE)
heatmap(normalized_counts[swapped_genes, ], scale = "row", col = colorRampPalette(c("blue", "white", "red"))(100))
