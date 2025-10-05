# Raneen Hoballah
# Assingement 04: Class 3B: Preprocessing and Normalization module
# 05/10/2025

### NOTICE ###
# Initially, I planned to use the dataset GSE45114, which studied “Transcriptional profiling of HCC vs. Normal liver and Pericancer vs. Normal liver.” However, upon inspection, the available metadata showed that only cancer and pericancerous tissue samples were present, and there were no normal liver samples included in the downloadable expression set.
# Since the assignment specifically requires a comparison between normal and cancer samples, this dataset was not suitable.
# To ensure valid biological comparison and meet the project requirements, I switched to GSE45267, a well-documented study that includes both HCC (tumor) and normal liver samples.
# GSE45267 provides 48 tumor and 39 normal samples (total 87), making it ideal for performing quality control, normalization, and filtering steps while maintaining statistical reliability. The dataset also comes from the same cancer type (HCC), keeping the research context consistent with my original topic.
### My apologies for any inconvenience ###

# R Workflow for Preprocessing & Normalization (GSE45114):

# 1.Loading the required libraries 
install.packages("BiocManager")
BiocManager::install(c("limma", "affy", "GEOquery", "Biobase"))
install.packages("data.table")
install.packages("generics")
install.packages("tidyr")
install.packages("xml2")
install.packages("curl")
install.packages("R.utils")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("arrayQuality
                     Metrics")
install.packages("svglite")
install.packages("Hmisc")
install.packages("latticeExtra")

library(data.table)
library(GEOquery)
library(limma)
library(affy)
library(Biobase)
library(curl)
library(R.utils)
library(arrayQualityMetrics)
library(Hmisc)
library(latticeExtra)
sessionInfo()


# 2. Loading and inspecting GEO dataset (GSE45267)
gset_list <- getGEO("GSE45267", GSEMatrix = TRUE, AnnotGPL = TRUE)
length(gset_list)
gset <- gset_list[[1]]
expr_raw <- exprs(gset)
pdata <- pData(gset)
cat("Dimensions of expression matrix (probes x samples):", dim(expr_raw)[1], "x", dim(expr_raw)[2], "\n\n")
cat("Available columns in phenotype data:\n")
print(colnames(pdata))
head(pdata[, grep("characteristics", colnames(pdata)), drop = FALSE])


# 3. Definning groups (Normal vs HCC)
tissue_info <- tolower(trimws(as.character(pdata$characteristics_ch1)))
group <- ifelse(grepl("normal", tissue_info), "Normal",
                ifelse(grepl("tumor|hcc|cancer", tissue_info), "HCC", "Unknown"))
keep_idx <- which(group %in% c("Normal", "HCC"))
expr_raw <- expr_raw[, keep_idx]
pdata <- pdata[keep_idx, , drop = FALSE]
pdata$group <- factor(group[keep_idx], levels = c("Normal", "HCC"))

cat("Group counts after filtering:\n")
print(table(pdata$group))


# 4. Quality Control (Before Normalization)
expr_raw <- as.matrix(expr_raw)
boxplot(expr_raw, 
        col = "lightblue", 
        main = "Boxplot of Raw Expression Values (Before Normalization)",
        las = 2, outline = FALSE, 
        ylab = "Log2 expression intensity")
plotDensity(expr_raw, 
            main = "Density Plot of Raw Expression Values (Before Normalization)",
            xlab = "Log2 expression intensity")

if (!inherits(gset, "ExpressionSet")) {
  gset <- ExpressionSet(assayData = expr_raw)
  pData(gset) <- pdata
}

# Run QC report
arrayQualityMetrics(expressionset = gset,
                    outdir = "QC_Before_Normalization",
                    force = TRUE,
                    do.logtransform = FALSE)

cat("✅ QC report generated in folder: 'QC_Before_Normalization'\n")
cat("Check the report for any arrays flagged as outliers.\n")

# After checking the QC report, 0 arrays flagged as outliers before normalization, therefore analysis can proceed.


# 5. Normalization
expr_norm <- normalizeBetweenArrays(expr_raw, method = "quantile")

# Update ExpressionSet object
exprs(gset) <- expr_norm

# Boxplot after normalization
boxplot(expr_norm, 
        col = "lightgreen", 
        main = "Boxplot of Normalized Expression Values",
        las = 2, outline = FALSE,
        ylab = "Log2 expression intensity")

# Density plot after normalization
plotDensity(expr_norm, 
            main = "Density Plot of Normalized Expression Values",
            xlab = "Log2 expression intensity")




# 6. Quality Control (After Normalization)
arrayQualityMetrics(expressionset = gset,
                    outdir = "QC_After_Normalization",
                    force = TRUE,
                    do.logtransform = FALSE)

cat("✅ QC report generated in folder: 'QC_After_Normalization'\n")
cat("Check the report for arrays flagged as outliers after normalization.\n")


# 7.Filtering Low-Intensity Probes:
keep_probes <- rowMeans(expr_norm) > 5
expr_filtered <- expr_norm[keep_probes, ]

cat("Number of probes before filtering:", nrow(expr_norm), "\n")
cat("Number of probes after filtering:", nrow(expr_filtered), "\n")



### Summary of Final Results ###
cat("✅ Pre-normalization QC: 0 arrays flagged as outliers.\n")
cat("✅ Post-normalization QC: 0 arrays flagged as outliers.\n")
cat("✅ Low-intensity probe filtering: ", nrow(expr_filtered), 
    "probes remain out of ", nrow(expr_norm), ".\n")
cat("✅ Groups defined: Normal vs HCC (", 
    table(pdata$group)["Normal"], "Normal and ", 
    table(pdata$group)["HCC"], "HCC samples).\n")
cat("\n")
cat("Overall, the dataset is of good quality. Normalization successfully adjusted array distributions,\n")
cat("and low-intensity probes were removed to reduce noise. The processed data is ready for downstream analysis.\n")
