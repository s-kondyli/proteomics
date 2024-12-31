# Load required libraries
library(dplyr)
library(tidyr)
library(VennDiagram)


#this script identifies common proteins between 2 samples having as input the pg.txt file from diann
# Convert columns that match the regex (basically the columns that contain LFQ normalized intensities) to numeric
columns_to_convert <- grep("^[t|s|r]", names(data), value = TRUE)
data[columns_to_convert] <- lapply(data[columns_to_convert], as.numeric)


# Function to clean up Gene Names from column "Genes" in DIANN 
data <- data %>%
  mutate(Genes = map(strsplit(Genes, ";"), 1))

#drop columns you don't need which are the first 3 columns and the 5th 
data <- data[, -c(1:3, 5)]

# Define a function to generate Venn diagram for two samples
generate_venn <- function(sample1, sample2, data) {
  # Get genes present in both samples and not 'Na'
  common_genes <- intersect(
    rownames(data[!is.na(data[, sample1]), ]),
    rownames(data[!is.na(data[, sample2]), ])
  )
  
  # Create a list of genes found in each sample
  genes_sample1 <- rownames(data[!is.na(data[, sample1]), ])
  genes_sample2 <- rownames(data[!is.na(data[, sample2]), ])
  
  # Generate Venn diagram
  venn_plot <- venn.diagram(
    x = list(Sample1 = genes_sample1, Sample2 = genes_sample2),
    category.names = c(sample1, sample2),
    filename = NULL
  )
  
  # Plot the Venn diagram
  grid.draw(venn_plot)
  
  # Return the common genes
  return(common_genes)
}

# Call the function with your desired samples -> adjust each time 
common_genes <- generate_venn("tip_10k_T_cells", "tip_1k_T_cells", data)

