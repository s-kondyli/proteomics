# This R script is designed for the LFQ analysis of mass spectrometry-based proteomics data from DIANN
# The script provides a systematic and reproducible analysis pipeline for proteomics data
# Lindeboom Lab - Cas Kranenburg feb 2024

## 1. **Environment Setup:**
# - Clears the environment and sets the seed for reproducibility.
# - Loads necessary libraries for data manipulation, visualization, and statistical analysis.

# Clean environment and set seed ----
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)
set.seed(123456)

# Libraries -----

library(DT)
library(MSnbase)
library(tibble)
library(SummarizedExperiment)
library(DEP)
library(shiny)
library(shinydashboard)
library(vsn)
library(limma)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(R.utils)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

# Before starting:

## 2.**Definition of parameters:**
# Requirement1: ExperimentalDesign.csv (colnames: "label" "condition" "replicate")
# Requirement2: raw data analysis output in tsv or txt format 


## 3. **Data Loading and Preprocessing:**
# - Reads proteomics output data along with experimental design files.
# - Merges the data and experimental design for downstream analysis.
# - Extracts gene names from protein descriptions, filters data, and identifies and resolves duplicate gene names.
  

# Load data -> Import the data using the "Import Dataset" from the environment
# Rename the raw data analysis PD output txt file to "data"
# Rename the  ExperimentalDesign.csv file to "ExperimentalDesign"



# Read all data 

ExperimentalDesign$replicate <- as.numeric(ExperimentalDesign$replicate )

#### temp code because of not all raw files were run -> delete from future experiments

# Find column names matching the regex pattern
columns_to_drop <- grep("^Abundances..Normalized...F(?:[4-6]|10|11)..Sample", names(data), value = TRUE)

# Drop the columns
data <- data[, !names(data) %in% columns_to_drop] # Now df contains all columns except those matching the regex pattern

#### end of temp code

# Convert columns that match the regex (basically the columns that contain maxLFQ intensities) to numeric
columns_to_convert <- grep("^Abundances..Normalized...F\\d..Sample", names(data), value = TRUE)
data[columns_to_convert] <- lapply(data[columns_to_convert], as.numeric)

# Function to extract Gene Name from Description
extract_gene_name <- function(description) {
  # Find the position of "GN=" in the Description
  gn_position <- regexpr("GN=", description)
  
  if (gn_position != -1) {
    # Extract the substring after "GN=" until the next space
    gene_name <- sub(".*GN=([^\\s]+).*", "\\1", description)
    gene_name <- sub(" PE=\\d+ SV=\\d+", "", gene_name)
    return(gene_name)
  } else {
    # If "GN=" is not found, return NA or handle accordingly
    return(NA)
  }
}

# Apply the function to the Description column and create a new column 'Gene.names'
data$Gene.names <- sapply(data$Description, extract_gene_name)

# Rename the column "Accession" to "Protein.IDs"
colnames(data)[colnames(data) == "Accession"] <- "Protein.IDs"


# filter the data if necessary (MQ)
# data <- filter(data, Reverse != "+", Potential.contaminant != "+")
dim(data)

#Verify any duplicate Gene.names
data$Genes.names %>% duplicated() %>% any()

# Make a table of duplicated gene names
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
# Verification to check duplicity. Aim: FALSE
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("^Abundances..Normalized...F\\d..Sample", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se(data_unique, LFQ_columns, ExperimentalDesign)
# Generate a SummarizedExperiment object by parsing condition information from the column names
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)



#Filter on missing values

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
# thr=1 <- Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
# thr=0 <- Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)


## 4. **Data Normalization and Imputation:**
# - Normalizes the data for downstream analysis..
# - Applies imputation techniques to handle missing values in the data

# Normalize the data
data_norm <- normalize_vsn(data_filt)
# meanSdPlot to verify the fitting of the normalization
# meanSdPlot(data_norm) # Can be skipped
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)



# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

#Imputation with algorithm

# Impute data for missing values- Imputation strategy
# Algorithm in fun should be one of:
# "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "man", "min", "zero", "mixed", "nbavg"
# data can be missing at random (MAR), for example if proteins are quantified in some replicates but not in others
# MAR data should be imputed with methods such as
# k-nearest neighbor (“knn”) or maximum likelihood (“MLE”) functions
data_imp <- impute(data_norm, fun = "knn")

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)


## 5. **Differential Enrichment Analysis:**
# - Performs differential enrichment analysis based on linear models and empirical Bayes statistics.
# - Identifies significant proteins based on user-defined criteria.

# Test every sample versus control
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1.5))

## 6. **Visualization of Results:**
# - Generates various plots, including PCA plots, correlation matrices, heatmaps, and volcano plots, to visualize the results of the differential enrichment analysis.

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

# Plot a volcano plot for the contrast -> change it each time according to your data
# in the plot_volcano function if you set adjusted = T it uses the adjusted p-value 
plot_volcano(dep, contrast = "plus.biotin.lysed.with.NP.40_vs_minus.biotin.lysed.with.NP.40", label_size = 2, add_names = TRUE, adjusted = T)  +
  labs(
  title = "Differential Enrichment plot of LS174T w4",
  subtitle = "Method: Astral DIA (30SPD) Ionoptics",
  x = (bquote(Log[2] ("fold-change"))),
  y = (bquote(-Log[10] ("adjusted p-value"))),
  caption = "LS174Tw4 cells +/- biotin lysed with NP-40 followed by strep pull-down"
) +
  #theme_classic() +
  theme (
    plot.title = element_text(size = 11),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 2, b = 1, l = 1))  # Adjust margin here
  ) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  

  
# Plot a barplot for a single protein -> change the protein name
plot_single(dep, proteins = c("PTK7", "CD47"))

# Plot a barplot for the single protein with the data centered
plot_single(dep, proteins = "PTK7", type = "centered")

# Plot a frequency plot of significant proteins for the different conditions -> figure out why it doesn't work
plot_cond(dep)

## 7. **Results Table:**
# - Produces a results table summarizing the outcomes of the differential enrichment analysis.
# - Calculates the number of significant proteins.

# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

colnames(data_results)

# Save the combined dataframe as a tsv txt file in the specified directory ->always change the file name 
data_folder <- "R:/Group Vermeulen/Lila/Mass_spec_results/V19_plus_minus_biotin"

write.table(data_results, file = paste0(data_folder, "/data_results_v19_plus_minus_biotin_PD.tsv"), sep = "\t", row.names = FALSE)



