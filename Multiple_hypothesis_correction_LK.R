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



## Multiple hypothesis correction using limma p.adjust function
# Input: data_results df from dep package

# Specify the regex pattern for columns to be dropped 
# because we change the correction some values may change the significance status
drop_pattern <- "*significant"  

# Find column names matching the drop pattern
columns_to_drop <- grep(drop_pattern, colnames(data_results), value = TRUE)

# Drop the columns from the dataframe
data_results <- data_results[, !colnames(data_results) %in% columns_to_drop]

####

# Specify the regex pattern for column names that contain the p values
regex_pattern <- "*p.val"  # Replace with your actual regex pattern

# Find column names matching the regex
matching_columns <- grep(regex_pattern, colnames(data_results), value = TRUE)

# Loop through matching columns and create new columns with adjusted p-values
for (col in matching_columns) {
  # Extract p-values from the original column
  p_values <- data_results[, col]
  
  # Apply p-value adjustment using limma's p.adjust function
  # You can choose other methods like "holm","hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  adjusted_p_values <- p.adjust(p_values, method = "bonferroni")  
  
  # Create a new column with adjusted p-values but keeping the rest of the name same and changing only the p.val for p.adj
  new_col_name <- sub("p.val", "p.adj", col)  # Replace "pval" with your actual column name pattern
  data_results[[new_col_name]] <- adjusted_p_values
}

### Save the data to be plotted 

# Save the combined dataframe as a tsv txt file in the specified directory ->always change the file name 
data_folder <- "R:/Group Vermeulen/Lila/Mass_spec_results/V18_plus_minus_biotin"

write.table(data_results, file = paste0(data_folder, "/GSEAdata_v18_plus_minus_biotin_DIANN.tsv"), sep = "\t", row.names = FALSE)