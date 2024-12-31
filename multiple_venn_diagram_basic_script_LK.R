
# Load required packages
library(dplyr)
library(VennDiagram)

# Up to 5 columns
# Choose the columns of interest ---> change each time
columns_of_interest <- c( "T_TILs_A_min_2", "T_TILs_A_min_3", "T_TILs_B_min_1", "T_TILs_B_min_2","T_TILs_B_min_3","T_TILs_C_min_1" )



# Function to clean up Gene Names from column "Genes" in DIANN 
data <- data %>%
  mutate(Genes = map(strsplit(Genes, ";"), 1))


# Convert columns that match the regex (basically the columns that contain maxLFQ intensities) to numeric
columns_to_convert <- grep("^[|T|c]", names(data), value = TRUE)
data[columns_to_convert] <- lapply(data[columns_to_convert], as.numeric)

# Binarize columns of interest
binarized_data <- data %>%
  mutate_at(vars(columns_of_interest), ~ if_else(is.na(.), 0, 1))


# Extract lists of proteins for each condition based on binarized columns
condition1_proteins <- binarized_data %>%
  filter( T_TILs_A_min_2 == 1) %>%
  pull(Genes)

condition2_proteins <- binarized_data %>%
  filter( T_TILs_A_min_3 == 1) %>%
  pull(Genes)

condition3_proteins <- binarized_data %>%
  filter( T_TILs_B_min_1 == 1) %>%
  pull(Genes)

condition4_proteins <- binarized_data %>%
  filter( T_TILs_B_min_2 == 1) %>%
  pull(Genes)

condition5_proteins <- binarized_data %>%
  filter( T_TILs_B_min_3 == 1) %>%
  pull(Genes)

condition6_proteins <- binarized_data %>%
  filter( T_TILs_C_min_1 == 1) %>%
  pull(Genes)

# Create a list of protein sets for the Venn diagram
protein_sets <- list(
  T_TILs_A_min_2 = condition1_proteins,
  T_TILs_A_min_3 = condition2_proteins,
  T_TILs_B_min_1 = condition3_proteins,
  T_TILs_B_min_2 = condition4_proteins,
  T_TILs_B_min_3 = condition5_proteins

)

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = protein_sets,
  category.names = c("A un 2", "A un 3", "B un 1","B un 2",  "B un 3"),
  filename = NULL,  # This prevents the diagram from being saved to a file
  output = TRUE
)

# Display the Venn diagram
grid.draw(venn.plot)
