# Raneen Hoballah
# Assignment 01: 02/08/2025

# Tasks: 
#  1. Set Working Directory: Done
# Create a new folder on your computer "AI_Omics_Internship_2025".
setwd("~/Downloads/AI_Omics_Internship_2025/Module_I")


# 2. Create Project Folder: Done 
# In RStudio, create a new project named "Module_I" in your "AI_Omics_Internship_2025" folder.
# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

dir.create("raw_data", showWarnings = FALSE)
dir.create("clean_data", showWarnings = FALSE)
dir.create("scripts", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)


# 3. Download "patient_info.csv" dataset from GitHub repository
# Downloading and loading the data:
url <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/main/patient_info.csv"
download.file(url, destfile = "raw_data/patient_info.csv", method = "curl")
patient_info <- read.csv ("raw_data/patient_info.csv")

# 4. Inspect the structure of the dataset using appropriate R functions
str(patient_info)
summary(patient_info)
head(patient_info)

# 5. Identify variables with incorrect or inconsistent data types.
# Based on the output of the previous code 
patient_info$age <- as.numeric(patient_info$age)
patient_info$bmi <- as.numeric(patient_info$bmi)
# Convert to factors
patient_info$gender <- as.factor(patient_info$gender)
patient_info$diagnosis <- as.factor(patient_info$diagnosis)
# Create a binary numeric variable for smoker, (additional for gender and diagnosis, just for the data to look nicer)
patient_info$smoker_binary <- ifelse(tolower(patient_info$smoker) == "yes", 1, 0)
patient_info$gender_binary <- ifelse(tolower(patient_info$gender) == "male", 1, 0)
patient_info$diagnosis_binary <- ifelse(tolower(patient_info$diagnosis) == "cancer", 1, 0)
# Deleting the extra columns (somker, gender and diagnosis)
patient_info$smoker <- NULL
patient_info$gender <- NULL
patient_info$diagnosis <- NULL
# Confirming the structure again
str(patient_info)

# 6. Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
write.csv(patient_info, "clean_data/patient_info_clean.csv", row.names = FALSE)

# 7. Save your R script in your script folder with name "class_Ib": Done
