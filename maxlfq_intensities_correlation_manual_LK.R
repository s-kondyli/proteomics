# Load necessary libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggcorrplot)
library(dplyr)

#Import data
proteomics_data = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all/all_filtered_2_unique_peptides/report.pg_peptides_filtered_all_PBMCs_DIANN_1.9.2 - Copy.tsv', stringsAsFactors = F, sep = '\t') # = pg.matrix

# Convert columns that match the regex (basically the columns that contain maxLFQ intensities) to numeric
columns_to_convert <- grep("wt|V|B", names(proteomics_data), value = TRUE)
proteomics_data[columns_to_convert] <- lapply(proteomics_data[columns_to_convert], as.numeric)

# Function to clean up Gene Names from column "Genes" in DIANN 
proteomics_data <- proteomics_data %>%
  mutate(Genes = map(strsplit(Genes, ";"), 1))

# Define a condition (replace "Condition1" with your desired condition)
condition_name <- "T"


# Filter experimental design for the condition
condition_samples <- ExperimentalDesign %>%
  filter(condition == condition_name) %>%
  pull(label)  # Extract sample labels for this condition

# Filter proteomics data for the condition's samples
proteomics_data_filtered <- proteomics_data %>%
  select(Genes, all_of(condition_samples))

# Remove rows with any NA values
proteomics_data_complete <- proteomics_data_filtered %>%
  filter(rowSums(is.na(.)) == 0)

# Apply log2 transformation
proteomics_data_log2 <- proteomics_data_complete %>%
  mutate(across(-Genes, ~ log2(.)))


# Step 1: Calculate the correlation matrix
correlation_matrix <- cor(
  proteomics_data_log2 %>% select(-Genes),
  use = "pairwise.complete.obs",
  method = "pearson"
)



# Generate the heatmap with more decimal precision
correlation_heatmap <- ggcorrplot(
  correlation_matrix,
  hc.order = TRUE,
  type = "full",
  lab = TRUE,                    # Add labels
  title = paste("Correlation Heatmap for Condition:", condition_name),
  colors = c("#ef8a62", "#f7f7f7", "#67a9cf"),
  digits = 3
)

# Plot the heatmap

print(correlation_heatmap)

