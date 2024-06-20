
# Load required packages
library(dplyr)
library(VennDiagram)

# Choose the columns of interest ---> change each time
columns_of_interest <- c("X1k.stimulated_vs_X1k.untreated_significant", "X2k.stimulated_vs_X2k.untreated_significant", 
                         "X5k.stimulated_vs_X5k.untreated_significant")

# Filter the dataframe to keep only the true values ----> skip this if you want a venn diagram for all proteins
filtered_data <- data %>%
  filter_at(vars(columns_of_interest), any_vars(. == TRUE))

# Extract lists of proteins for each condition -----> change each time according to your conditions
condition1_proteins <- filtered_data %>%
  filter(X1k.stimulated_vs_X1k.untreated_significant == TRUE) %>%
  pull(name) #name is the column that contains the proteins names

condition2_proteins <- filtered_data %>%
  filter(X2k.stimulated_vs_X2k.untreated_significant == TRUE) %>%
  pull(name)

condition3_proteins <- filtered_data %>%
  filter(X5k.stimulated_vs_X5k.untreated_significant == TRUE) %>%
  pull(name) 

# Create a list of protein sets for the Venn diagram
protein_sets <- list(
  X1k.stimulated_vs_X1k.untreated_significant = condition1_proteins,
  X2k.stimulated_vs_X2k.untreated_significant = condition2_proteins,
  X5k.stimulated_vs_X5k.untreated_significant = condition3_proteins
)

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = protein_sets,
  category.names = c("1k T cells", "2k T cells", "5K T cells"),
  filename = NULL,  # This prevents the diagram from being saved to a file
  output = TRUE
)

# Display the Venn diagram
grid.draw(venn.plot)






# Function to count significant proteins in a specific condition
count_significant_proteins <- function(data, condition_column) {
  count <- data %>%
    filter(!!sym(condition_column) == TRUE) %>%
    nrow()
  return(count)
}

# Example usage: count significant proteins for Condition1
condition_column <- "X5k.stimulated_vs_X5k.untreated_significant"
significant_proteins_count <- count_significant_proteins(data, condition_column)

# Display the count of significant proteins for the specified condition
significant_proteins_count
