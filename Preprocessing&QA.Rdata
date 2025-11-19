# Group Project
# Raneen Hoballah: Preprocessing and QA


# Install required packages (if not installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "affy", "simpleaffy", "arrayQualityMetrics", "limma"))

library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(limma)


# Dowenload and extract the data (choose the directory that you like)
setwd("~/Downloads")
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15852/suppl/GSE15852_RAW.tar"
download.file(url, destfile = "~/Downloads/GSE15852_RAW.tar", mode = "wb", method = "curl")

# Unzip the files
untar("GSE15852/GSE15852_RAW.tar", exdir = "GSE15852/CEL")
cel_files <- list.files("GSE15852/CEL", pattern = "gz$", full.names = TRUE)
sapply(cel_files, R.utils::gunzip, overwrite = TRUE)
file.info("~/Downloads/GSE15852_RAW.tar")

# Chceking if everything is alright
setwd("~/Downloads")
dir.create("GSE15852_CEL", showWarnings = FALSE)
untar("GSE15852_RAW.tar", exdir = "GSE15852_CEL")
list.files("GSE15852_CEL")[1:10]

# Read and preprocess the files
raw_data <- ReadAffy(celfile.path = "GSE15852_CEL")
norm_data <- rma(raw_data) # RMA normalization
exprs_mat <- exprs(norm_data) # Extract expression matrix
dim(exprs_mat)

# Box Plot
boxplot(exprs_mat, main = "RMA-normalized expression", las = 2, col = "lightblue")

# Box Plot Interpretation:
# Each box represents one CEL file (sample). The median (black line) and interquartile range (IQR) of expression values are plotted.
# All boxes are nicely aligned around the same median (~10). The spread (IQR) is consistent across samples. There are no extreme outliers or shifted distributions. This indicates successful normalization — the data are now comparable across samples.

# Plot the Density of all samples
boxplot(exprs_mat, main = "RMA-normalized expression", las = 2, col = "lightblue")
plotDensities(exprs_mat, main = "Density plot of RMA-normalized data")
arrayQualityMetrics(norm_data, outdir = "QA_Report_GSE15852", force = TRUE)

#### Everything in the reports looks good, this is the end of QA ####

# Assign Tumor vs Normal Labels
gse <- getGEO("GSE15852", GSEMatrix = TRUE)[[1]]

# Extract phenotype data
pheno <- pData(gse)
head(pheno$characteristics_ch1)

# Create normal/tumor factor
group <- ifelse(grepl("tumor", pheno$characteristics_ch1, ignore.case = TRUE),
                "Tumor", "Normal")
group <- factor(group)
table(group)

# Chacking and Verifying 
cel_names <- gsub(".CEL$", "", list.celfiles("GSE15852_CEL", full.names = FALSE))
head(cel_names)
head(rownames(pheno)) # Compare with GEO sample names


# PCA on normalized expression
pca <- prcomp(t(exprs_mat), scale. = TRUE)
plot(pca$x[,1], pca$x[,2],
     col = ifelse(group == "Tumor", "red", "blue"),
     pch = 19,
     main = "PCA of GSE15852 samples",
     xlab = "PC1", ylab = "PC2")
legend("topright", legend = levels(group), col = c("blue","red"), pch = 19)





