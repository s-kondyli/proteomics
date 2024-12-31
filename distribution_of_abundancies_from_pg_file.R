# Load required libraries
library(ggplot2)
library(reshape2)

# Read the data
data <- read.delim("R:/Group Vermeulen/Lila/Mass_spec_results/V65_LK_Astral/report.pg_peptides_filtered - Copy.tsv", sep = "\t", header = TRUE)



# Select only the columns with abundance values (assuming they start after the first few columns)
abundance_data <- data[, grep("GFP|colo", colnames(data))]

# Convert all abundance columns to numeric, ensuring non-numeric values are handled
abundance_data <- as.data.frame(lapply(abundance_data, function(x) {
  # Convert to numeric after removing commas or special characters
  as.numeric(gsub(",", "", as.character(x)))
}))

# Apply log2 transformation to abundance data, ignoring zero or negative values to avoid log(0) issues
abundance_data <- as.data.frame(lapply(abundance_data, function(x) log2(x[x > 0])))

# Melt the data to a long format for easier plotting
melted_data <- melt(abundance_data, variable.name = "Sample", value.name = "Log2_Abundance")

# Remove any NA values from the data
melted_data <- melted_data[!is.na(melted_data$Log2_Abundance), ]





# Plot the distribution of log2-transformed abundances
ggplot(melted_data, aes(x = Log2_Abundance)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ Sample, scales = "free_x") +
  scale_x_continuous(breaks = seq(5, 25, by = 3)) +  # Custom breaks for x-axis
  theme_minimal() +
  labs(title = "Distribution of Log2-Transformed Intensities Across Samples",
       x = "Log2 Intensity",
       y = "Count") 







#----------------------------------- if you want to annotate a certain protein
# Define the list of proteins of interest based on 'Protein.Names' column
proteins_of_interest <- c("SMAD4_HUMAN")  # Replace with actual protein names of interest

# Select only the columns with abundance values (assuming they start after the first few columns)
abundance_data <- data[, grep("GFP|colo", colnames(data))]

# Convert all abundance columns to numeric, ensuring non-numeric values are handled
abundance_data <- as.data.frame(lapply(abundance_data, function(x) {
  # Convert to numeric after removing commas or special characters
  as.numeric(gsub(",", "", as.character(x)))
}))

# Apply log2 transformation to abundance data, ignoring zero or negative values to avoid log(0) issues
abundance_data <- as.data.frame(lapply(abundance_data, function(x) log2(x[x > 0])))

# Melt the data to a long format for easier plotting
melted_data <- melt(abundance_data, variable.name = "Sample", value.name = "Log2_Abundance")

# Remove any NA values from the data
melted_data <- melted_data[!is.na(melted_data$Log2_Abundance), ]

# Extract  abundance values for proteins of interest
protein_interest_data <- data[data$Protein.Names %in% proteins_of_interest, ]

# Extract abundance columns for proteins of interest and log2-transform them
protein_interest_abundance <- protein_interest_data[, grep("GFP|colo", colnames(protein_interest_data))]
protein_interest_abundance_log2 <- as.data.frame(lapply(protein_interest_abundance, function(x) log2(as.numeric(gsub(",", "", as.character(x))))))

# Melt protein interest data for easier plotting
protein_interest_data_melted <- melt(protein_interest_abundance_log2, variable.name = "Sample", value.name = "Log2_Abundance")



# Plot the distribution of log2-transformed abundances with protein of interest overlay
ggplot(melted_data, aes(x = Log2_Abundance)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ Sample, scales = "free_x") +
  geom_point(data = protein_interest_data_melted, aes(x = Log2_Abundance, y = 0), 
             color = "red", size = 2, shape = 17) +  # Overlay protein of interest
  scale_x_continuous(breaks = seq(5, 25, by = 3)) +  # Custom breaks for x-axis
  theme_minimal() +
  labs(title = "Distribution of Log2-Transformed Intensities Across Samples",
       x = "Log2 Intensity",
       y = "Count") +
  theme(legend.position = "bottom")




