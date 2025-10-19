# Raneen Hoballah
# Assignment 05: Class 3CB: Differential Gene Expression Analysis
# 19/10/2025

# Dataset: GSE45267

# 1. Install & load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("GEOquery","limma","AnnotationDbi","hgu133plus2.db","Biobase")
pkgs_cran <- c("pheatmap","ggplot2","dplyr","tibble","matrixStats")

for (p in c(pkgs_bioc, pkgs_cran)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% pkgs_bioc) BiocManager::install(p, ask = FALSE, update = FALSE) else install.packages(p)
  }
  library(p, character.only = TRUE)
}


# 2. Create results folder
if (!dir.exists("Results")) dir.create("Results")


# 3. Download GSE45267 from GEO
gset_list <- getGEO("GSE45267", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset_list[[1]]    # pick the first ExpressionSet


# 4. Extract expression matrix & phenotype
expr_raw <- exprs(gset)   
pdata <- pData(gset)
cat("Expression matrix dimensions (probes x samples):", dim(expr_raw), "\n")
cat("Phenotype columns available:\n"); print(colnames(pdata))


# 5. Define groups: Normal vs HCC
char_cols <- grep("characteristics|title|source|tissue", tolower(colnames(pdata)), value = TRUE)
found <- NULL
for (cc in char_cols) {
  v <- tolower(as.character(pdata[[cc]]))
  if (any(grepl("normal|tumor|hcc|cancer", v))) { found <- cc; break }
}
if (is.null(found)) {
  
  found <- if ("characteristics_ch1" %in% colnames(pdata)) "characteristics_ch1" else colnames(pdata)[1]
}
cat("Using phenotype column for grouping:", found, "\n")
tissue_info <- tolower(trimws(as.character(pdata[[found]])))
group <- ifelse(grepl("normal", tissue_info), "Normal",
                ifelse(grepl("tumor|hcc|cancer|carcinoma|primary tumour|primary tumor", tissue_info), "HCC", "Other"))
keep_idx <- which(group %in% c("Normal","HCC"))
expr_raw <- expr_raw[, keep_idx]
pdata <- pdata[keep_idx, , drop = FALSE]
pdata$group <- factor(group[keep_idx], levels = c("Normal","HCC"))
cat("Group counts:\n"); print(table(pdata$group))

# 6. Normalization
expr_norm <- normalizeBetweenArrays(as.matrix(expr_raw), method = "quantile")


# 7. Filter low-intensity probes
keep_probes <- rowMeans(expr_norm) > 5  
expr_filtered <- expr_norm[keep_probes, ]
cat("Probes before filtering:", nrow(expr_norm), "after filtering:", nrow(expr_filtered), "\n")


# 8. Map probe IDs to gene symbols using hgu133plus2.db
library(hgu133plus2.db)
probe_ids <- rownames(expr_filtered)
symbols <- AnnotationDbi::mapIds(hgu133plus2.db,
                                 keys = probe_ids,
                                 column = "SYMBOL",
                                 keytype = "PROBEID",
                                 multiVals = "first")  
map_df <- data.frame(PROBEID = probe_ids,
                     SYMBOL = as.character(symbols),
                     stringsAsFactors = FALSE)
mapped_non_na <- map_df[!is.na(map_df$SYMBOL) & map_df$SYMBOL != "", ]
probe_count_per_gene <- table(mapped_non_na$SYMBOL)
genes_with_multiple_probes <- names(probe_count_per_gene[probe_count_per_gene > 1])
n_genes_multi <- length(genes_with_multiple_probes)
n_probes_mapped <- nrow(mapped_non_na)
n_probes_total <- length(probe_ids)
cat("Probes with any gene symbol mapped:", n_probes_mapped, "out of", n_probes_total, "\n")
cat("Number of genes with multiple probes mapped:", n_genes_multi, "\n")

# 9. Handle duplicates: collapse probes by gene symbol
# Strategy: for each gene, I'll keep the probe with the highest mean expression across samples
expr_mapped <- expr_filtered[rownames(expr_filtered) %in% mapped_non_na$PROBEID, ]
probe_to_symbol <- mapped_non_na$SYMBOL
names(probe_to_symbol) <- mapped_non_na$PROBEID
probe_means <- rowMeans(expr_mapped)
probe_info <- data.frame(PROBEID = rownames(expr_mapped), SYMBOL = probe_to_symbol[rownames(expr_mapped)], meanExpr = probe_means, stringsAsFactors = FALSE)

library(dplyr)
rep_probes <- probe_info %>%
  group_by(SYMBOL) %>%
  arrange(desc(meanExpr)) %>%
  slice(1) %>%
  ungroup()
expr_collapsed <- expr_mapped[rep_probes$PROBEID, ]
rownames(expr_collapsed) <- rep_probes$SYMBOL
cat("Collapsed expression matrix dimensions (genes x samples):", dim(expr_collapsed), "\n")

# 10. Differential expression with limma: HCC vs Normal
group <- pdata$group
design <- model.matrix(~0 + group)
colnames(design) <- c("Normal","HCC")
contrast <- makeContrasts(HCC_vs_Normal = HCC - Normal, levels = design)
fit <- lmFit(expr_collapsed, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
res_all <- topTable(fit2, coef = "HCC_vs_Normal", number = Inf, adjust.method = "BH", sort.by = "P")

# 11. Set DEG thresholds
adj_p_cut <- 0.05
logfc_cut <- 1
res_all$Gene <- rownames(res_all)
write.csv(res_all, file = "Results/all_DEGs_full_table.csv", row.names = FALSE)
res_up <- res_all %>% filter(adj.P.Val <= adj_p_cut & logFC >= logfc_cut)
res_down <- res_all %>% filter(adj.P.Val <= adj_p_cut & logFC <= -logfc_cut)
write.csv(res_all, file = "Results/DEG_full_table.csv", row.names = FALSE)
write.csv(res_up, file = "Results/upregulated_DEGs.csv", row.names = FALSE)
write.csv(res_down, file = "Results/downregulated_DEGs.csv", row.names = FALSE)
cat("DEG files written to Results/ (all, upregulated, downregulated)\n")
cat("Upregulated genes:", nrow(res_up), " Downregulated genes:", nrow(res_down), "\n")


### PLOTS ###
# 12. Volcano plot (ggplot2)
library(ggplot2)
res_all$label <- ifelse(res_all$adj.P.Val <= adj_p_cut & abs(res_all$logFC) >= logfc_cut,
                        ifelse(res_all$logFC > 0, "Up","Down"), "NotSig")
res_all$negLogP <- -log10(res_all$P.Value)

p <- ggplot(res_all, aes(x = logFC, y = negLogP, color = label)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey")) +
  geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(adj_p_cut), linetype = "dashed") +
  labs(title = "Volcano plot: HCC vs Normal", x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()

ggsave("Results/volcano.png", plot = p, width = 7, height = 6, dpi = 300)
cat("Volcano plot saved to Results/volcano.png\n")

# 13. Heatmap of top 25 DEGs
topN <- 25
top_genes <- head(res_all$Gene[!is.na(res_all$adj.P.Val)], topN)
expr_top <- expr_collapsed[top_genes, , drop = FALSE]
scale_rows <- function(x) {
  m <- rowMeans(x, na.rm = TRUE)
  s <- matrixStats::rowSds(as.matrix(x), na.rm = TRUE)
  s[s == 0] <- 1
  (x - m) / s
}
expr_top_z <- scale_rows(expr_top)
annotation_col <- data.frame(Group = pdata$group)
rownames(annotation_col) <- colnames(expr_top_z)

pheatmap::pheatmap(expr_top_z,
                   filename = "Results/heatmap_top25.png",
                   main = "Top 25 DEGs (HCC vs Normal)",
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   annotation_col = annotation_col,
                   show_rownames = TRUE, show_colnames = FALSE,
                   width = 8, height = 10)

cat("Heatmap saved to Results/heatmap_top25.png\n")

# 14. Write a short summary file
summary_lines <- c(
  "Class 3CB: Differential Expression Analysis - Summary",
  paste0("Dataset: GSE45267 (GPL570, Affymetrix Human Genome U133 Plus 2.0)"),
  "",
  paste0("Total probes after initial filtering: ", n_probes_total),
  paste0("Probes mapped to gene symbols: ", n_probes_mapped),
  paste0("Number of genes represented after collapsing duplicates: ", nrow(expr_collapsed)),
  paste0("Number of genes that had multiple probes mapped: ", n_genes_multi),
  "Duplicate probe handling strategy: For genes with multiple probes, the probe with the highest mean expression across all samples was chosen as the representative probe.",
  "",
  paste0("DE contrast performed: HCC_vs_Normal (HCC minus Normal)"),
  paste0("DEG thresholds used: adjusted p-value <= ", adj_p_cut, ", |log2FC| >= ", logfc_cut),
  paste0("Total genes tested (after collapse): ", nrow(res_all)),
  paste0("Number of upregulated genes (adj.P.Val <= ", adj_p_cut, " & logFC >= ", logfc_cut, "): ", nrow(res_up)),
  paste0("Number of downregulated genes (adj.P.Val <= ", adj_p_cut, " & logFC <= -", logfc_cut, "): ", nrow(res_down)),
  "",
  "Files generated:",
  "- Results/DEG_full_table.csv (complete limma results)",
  "- Results/upregulated_DEGs.csv",
  "- Results/downregulated_DEGs.csv",
  "- Results/volcano.png",
  "- Results/heatmap_top25.png",
  "",
  paste0("Script run date: ", Sys.Date())
)

writeLines(summary_lines, con = "Results/summary.txt")
cat("Summary written to Results/summary.txt\n")

# 14. Last Step: End
cat("All steps complete. Check the Results/ folder for CSVs, PNGs, and summary.txt\n")

