# Load necessary libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggcorrplot)

#Creating a correlation plot to compare protein abundances among replicates of 1 condition
#is useful to assess the reproducibility of your mass spectrometry data
# Input: ExperimentalDesign & pg.matrix 

#Import data
proteomics_data = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all/report.pg_matrix - Copy.tsv', stringsAsFactors = F, sep = '\t') # = pg.matrix

# Convert columns that match the regex (basically the columns that contain maxLFQ intensities) to numeric
columns_to_convert <- grep("wt|V|B", names(proteomics_data), value = TRUE)
proteomics_data[columns_to_convert] <- lapply(proteomics_data[columns_to_convert], as.numeric)


# Extract only sample columns (assume sample columns start after metadata)
sample_columns <- grep("V|B", colnames(proteomics_data), value = TRUE)

# Filter out excluded samples based on experimental design
excluded_samples <- setdiff(sample_columns, ExperimentalDesign$label)
proteomics_data_filtered <- proteomics_data %>%
  select(-one_of(excluded_samples))

# Merge design information to create annotations
annotated_data <- ExperimentalDesign %>%
  filter(label %in% colnames(proteomics_data_filtered)) %>%
  mutate(label = factor(label, levels = colnames(proteomics_data_filtered)))

# Use only the remaining sample columns after exclusion
remaining_sample_columns <- colnames(proteomics_data_filtered)[grepl("^V", colnames(proteomics_data_filtered))]

# Reshape proteomics data to long format
long_data <- proteomics_data_filtered %>%
  select(Protein.Group, all_of(remaining_sample_columns)) %>%
  melt(id.vars = "Protein.Group", variable.name = "label", value.name = "abundance")

# Merge annotation data with long format proteomics data
annotated_long_data <- long_data %>%
  inner_join(annotated_data, by = "label")




# Function to create heatmap for each condition with log2 transformation
create_heatmap <- function(condition, data) {
  # Filter data for the given condition
  condition_data <- data %>% filter(condition == !!condition)
  
  # Remove proteins with NA values in any replicate
  complete_proteins <- condition_data %>%
    group_by(Protein.Group) %>%
    filter(all(!is.na(abundance))) %>%
    ungroup()
  
  # Check if Protein.Group exists
  if (!"Protein.Group" %in% colnames(complete_proteins)) {
    stop("Error: 'Protein.Group' column is missing in the filtered data.")
  }
  
  # Reshape to wide format
  wide_data <- complete_proteins %>%
    select(Protein.Group, label, abundance) %>%
    pivot_wider(names_from = label, values_from = abundance)
  
  # Apply log2 transformation 
  wide_data <- wide_data %>%
    mutate(across(-Protein.Group, ~ log2(.)))  # Apply log2 transformation to abundance values
  
  # Remove Protein.Group for correlation calculation
  wide_data_matrix <- wide_data %>% select(-Protein.Group)
  
  # Calculate correlation matrix
  correlation_matrix <- cor(wide_data_matrix, use = "pairwise.complete.obs", method = "pearson")
  
  # Generate heatmap
  heatmap <- ggcorrplot(
    correlation_matrix,
    hc.order = TRUE, # Hierarchical clustering
    type = "full",
    lab = TRUE,
    title = paste("Correlation Heatmap for Condition:", condition),
    colors = c("red", "white", "blue")
  )
  
  return(heatmap)
}



# Apply the function to each condition
unique_conditions <- unique(annotated_data$condition)
heatmaps <- lapply(unique_conditions, create_heatmap, data = annotated_long_data)

print(heatmaps)
