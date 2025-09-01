# Raneen Hoballah
# Assingment 02:  31/08/2025

# Q1: Write the function classify_gene()
classify_gene <- function(logFC, padj) {
  if (!is.numeric(logFC) | !is.numeric(padj)) {
    return("Not_Significant")  
  } else if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}


# Q2: Load the data
# Load datasets
DEGs_data_1 <- read.csv("/Users/raneenhob/Downloads/AI_Omics_Internship_2025/Module_I/raw_data/DEGs_Data_1.csv", stringsAsFactors = FALSE)
DEGs_data_2 <- read.csv("/Users/raneenhob/Downloads/AI_Omics_Internship_2025/Module_I/raw_data/DEGs_Data_2.csv", stringsAsFactors = FALSE)
# Checking the first few rows just to confirm
head(DEGs_data_1)
head(DEGs_data_2)


# Q3: Handle missing values
# Checking if there is missing data: there is 2 ways to do so
sum(is.na(DEGs_data_1$padj))
sum(is.na(DEGs_data_2$padj))
colSums(is.na(DEGs_data_1))
colSums(is.na(DEGs_data_2))
# Based on the output, there is no missing data therefore no need to do any furthere changes


# Q4: Apply the classify_gene function row by row
DEGs_data_1$status <- NA
for (i in 1:nrow(DEGs_data_1)) {
  DEGs_data_1$status[i] <- classify_gene(DEGs_data_1$logFC[i], DEGs_data_1$padj[i])
}
DEGs_data_2$status <- NA
for (i in 1:nrow(DEGs_data_2)) {
  DEGs_data_2$status[i] <- classify_gene(DEGs_data_2$logFC[i], DEGs_data_2$padj[i])
}


# Q5: Save the processed datasets
write.csv(DEGs_data_1, "/Users/raneenhob/Downloads/AI_Omics_Internship_2025/Module_I/results/DEGs_data_1_processed.csv", row.names = FALSE)
write.csv(DEGs_data_2, "/Users/raneenhob/Downloads/AI_Omics_Internship_2025/Module_I/results/DEGs_data_2_processed.csv", row.names = FALSE)


# Q6: Summarize results
cat("DEGs_data_1 counts:\n")
print(table(DEGs_data_1$status))
cat("\nDEGs_data_2 counts:\n")
print(table(DEGs_data_2$status))
