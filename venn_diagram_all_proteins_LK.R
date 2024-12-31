# Load required libraries
library(VennDiagram)
library(tidyverse)

# Step 1: Load the data
data = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/LK_V48_Paulien/report.pg_matrix - Copy.tsv', stringsAsFactors = F, sep = '\t') # = pg.matrix

# Step 2: Process the experimental design
colnames(ExperimentalDesign) <- c("label", "condition", "replicate")

# Step 3: Group replicates by condition
condition <- "tip" # Change this to the condition of interest
replicates <- ExperimentalDesign %>%
  dplyr::filter(condition == !!condition) %>% # Use dplyr::filter explicitly
  pull(label)

# Step 4: Extract gene lists for each replicate
genes_by_replicate <- lapply(replicates, function(replicate) {
  data %>%
    dplyr::filter(!is.na(!!sym(replicate))) %>% # Use dplyr::filter explicitly
    pull(Genes)
})
names(genes_by_replicate) <- replicates

# Step 5: Create the Venn Diagram
venn.plot <- venn.diagram(
  x = genes_by_replicate,
  category.names = names(genes_by_replicate),
  filename = NULL, # Set to NULL to render to R graphics device
  output = TRUE,
  col = "black",
  fill = c("#1b9e77", "#d95f02", "#7570b3", "#a6611a"), # Adjust number of colors depending on the numbers of samples #a6611a #7570b3
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  cat.col = "black",
  cat.cex = 1.2,
  main ="Overlap of proteins identified in all samples sorted in tips"
)

# Plot the Venn Diagram
grid.newpage()
grid.draw(venn.plot)

