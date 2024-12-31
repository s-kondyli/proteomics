# Load necessary libraries
library(tidyverse)
library(VennDiagram)



# Define contrasts
contrasts <- c("V88_T_vs_V88_B", "V85_T_vs_V85_B", "V93_T_vs_V93_B")

#Ratio threshold for significance
ratio_threshold <- log2(1.5)


#Define significant proteins
significant_proteins <- list()

for (contrast in contrasts) {
  significant_proteins[[contrast]] <- data %>% 
    filter(
      !!sym(paste0(contrast, "_p.adj")) < 0.05 &
        (!!sym(paste0(contrast, "_ratio")) > ratio_threshold | 
           !!sym(paste0(contrast, "_ratio")) < -ratio_threshold)
    ) %>% 
    pull(name)
}

# Create Venn diagram for significant proteins
venn_significant <- venn.diagram(
  x = significant_proteins,
  category.names = contrasts,
  filename = NULL,
  fill = c("#1b9e77", "#d95f02", "#7570b3"),
  alpha = 0.5,
  main = "Significant Proteins (p.adj < 0.05 & ratio > |log2(1.5)|)"
)

# Plot the Venn diagram for significant proteins
grid.newpage()
grid.draw(venn_significant)

