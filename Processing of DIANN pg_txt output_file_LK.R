## 1. **Environment Setup:**
# - Clears the environment and sets the seed for reproducibility.
# - Loads necessary libraries for data manipulation, visualization, and statistical analysis.

# Clean environment and set seed ----
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)
set.seed(123456)

# Libraries -----
# Install and load the required packages
#install.packages(c("tidyverse", "tibble", "RColorBrewer", "ggrepel", "SummarizedExperiment", "shiny", "shinydashboard", "R.utils"))
# BiocManager::install(c("DT", "MSnbase", "DEP", "limma", "vsn"))
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


## 2.**Definition of parameters:**
# Requirement1: ExperimentalDesign.csv (colnames: "label" "condition" "replicate")
# Requirement2: raw data analysis output in tsv or txt format 


## 3. **Data Loading and Preprocessing:**
# - Reads proteomics output data along with experimental design files.
# - Merges the data and experimental design for downstream analysis.
# - Extracts gene names from protein descriptions, filters data, and identifies and resolves duplicate gene names.

# Load data 
# Import the data using the "Import Dataset" from the environment
# Rename the raw data analysis output file to "data"
# Rename the  ExperimentalDesign.csv file to "ExperimentalDesign"



# Read all data 

ExperimentalDesign$replicate <- as.numeric(ExperimentalDesign$replicate )

# Convert columns that match the regex (basically the columns that contain maxLFQ intensities) to numeric
columns_to_convert <- grep("^[C|u]", names(data), value = TRUE)
data[columns_to_convert] <- lapply(data[columns_to_convert], as.numeric)


# Function to clean up Gene Names from column "Genes" in DIANN 
data <- data %>%
  mutate(Genes = map(strsplit(Genes, ";"), 1))

# filter the data if necessary (MQ)
# data <- filter(data, Reverse != "+", Potential.contaminant != "+")
dim(data)

#Verify any duplicate Gene.names
data$Genes %>% duplicated() %>% any()

# Make a table of duplicated gene names
data %>% group_by(Genes) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
# Verification to check duplicity. Aim: FALSE
data_unique <- make_unique(data, "Genes", "Protein.Ids", delim = ";")

# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("^[C|u]", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se(data_unique, LFQ_columns, ExperimentalDesign)
# Generate a SummarizedExperiment object by parsing condition information from the column names
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)


# Filter for proteins that are identified in all replicates of at least one condition
# thr=1 <- Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
# thr=0 <- Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

#Imputation with algorithm

# Impute data for missing values- Imputation strategy
# Algorithm in fun should be one of:
# "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "man", "min", "zero", "mixed", "nbavg"
# data can be missing at random (MAR), for example if proteins are quantified in some replicates but not in others
# MNAR (missing not at random ->relevant to use) data should be imputed with left-censored methods such as MinProb, min, zero
data_imp <- impute(data_norm, fun = "MinProb")


# Test every sample versus control
data_diff_all_contrasts <- test_diff(data_imp, type = "all")



# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1.5))


## 7. **Results Table:**
# - Produces a results table summarizing the outcomes of the differential enrichment analysis.
# - Calculates the number of significant proteins.

# Generate a results table
data_results <- get_results(dep)

#Filter out contaminants e.g. keratins etc

#Create a list of contaminants form the common_contaminants txt file
# Read the file line by line
lines <- readLines("R:/Group Vermeulen/Lila/Mass_spec_results/common_contaminants.txt")

# Initialize an empty list to store the results
sequences <- list()

# Loop through each line
for (line in lines) {
  # Check if the line starts with ">"
  if (startsWith(line, ">")) {
    # Extract the part after ">" and before the next space
    sequence_id <- gsub("^>\\s*([^ ]+).*", "\\1", line)
    # Append the extracted part to the list
    sequences <- c(sequences, list(sequence_id))
  }
}
# Remove rows based on sequences list
data_results <- data_results[!(data_results$ID %in% sequences), ]


# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

colnames(data_results)

# Save the combined dataframe as a tsv txt file in the specified directory ->always change the file name 
data_folder <- "R:/Group Vermeulen/Lila/Mass_spec_results/V26_Paulien"

write.table(data_results, file = paste0(data_folder, "/data_results_V26_DIANN_.tsv"), sep = "\t", row.names = FALSE)



### Visualization 

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)



# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

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
plot_volcano(dep, contrast = "wt_vs_PDL2_ko", label_size = 2, add_names = TRUE, adjusted = T) +
  labs(
    title = "Differential Enrichment plot",
    subtitle = "Method: Astral DIA (30SPD) Ionoptics",
    x = (bquote(Log[2] ("fold-change"))),
    y = (bquote(-Log[10] ("adjusted p-value"))),
    caption = "."
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
plot_single(dep, proteins = "TFRC", type = "centered")

# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)
