# Load required libraries
library(dplyr)
library(ggplot2)


# Keep only the relevant columns with the iBAQ values
iBAQ_info <- iBAQ_info %>%
  select(Protein.Group, Genes, contains("_iBAQ"))

# Define regex pattern and identify the column
sample_column <- grep("CUPS_4", colnames(iBAQ_info), value = TRUE)
if (length(sample_column) != 1) stop("The regex matched zero or multiple columns.")

# Filter out rows with infinite iBAQ values before doing anything
iBAQ_info <- iBAQ_info %>%
  filter(!is.infinite(.data[[sample_column]]))  # Exclude rows with infinite iBAQ values

# Filter, log2 transform, and exclude infinite values
top_100_data <- iBAQ_info %>%
  select(Genes, all_of(sample_column)) %>% # Keep only Genes and the selected sample column
  filter(!is.na(.data[[sample_column]])) %>% # Remove NA values
  mutate(log2_iBAQ = log2(.data[[sample_column]] + 1)) %>% # Log2 transform
  arrange(desc(log2_iBAQ)) %>% # Sort in descending order
  slice(1:100) # Select top 100 proteins


#Create the bar plot
ggplot(top_100_data, aes(x = reorder(Genes, -log2_iBAQ), y = log2_iBAQ)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    title = paste("Top 100 Most Abundant Proteins in CUP 4"),
    x = "Genes",
    y = "Log2 Transformed iBAQ Values"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels