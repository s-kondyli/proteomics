# Load required libraries
library(dplyr)
library(tidyr)
library(VennDiagram)
library(DT)
library(MSnbase)
library(tibble)
library(shiny)
library(shinydashboard)
library(vsn)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(R.utils)


# Convert columns that match the regex (basically the columns that contain maxLFQ intensities) to numeric
columns_to_convert <- grep("^[|t|s|r]", names(data), value = TRUE)
data[columns_to_convert] <- lapply(data[columns_to_convert], as.numeric)


# Function to clean up Gene Names from column "Genes" in DIANN 
data <- data %>%
  mutate(Genes = map(strsplit(Genes, ";"), 1))

# Function to clean up Gene Names from column "Protein.Ids" in DIANN 
data <- data %>%
  mutate(Protein.Ids = map(strsplit(Protein.Ids, ";"), 1))

#Filter out contaminants e.g. keratins etc

#Create a list of contaminants form the common_contaminants txt file
# Read the file line by line
lines <- readLines("R:/Group Vermeulen/Lila/Mass_spec_results/common_contaminants.txt")

# Initialize an empty list to store the results
sequences <- list()

# Loop through each line
for (line in lines) {
  # Check if the line starts with ">"
  if (startsWith(line, ">")) {
    # Extract the part after ">" and before the next space
    sequence_id <- gsub("^>\\s*([^ ]+).*", "\\1", line)
    # Append the extracted part to the list
    sequences <- c(sequences, list(sequence_id))
  }
}
# Remove rows based on sequences list
data <- data[!(data$Protein.Ids %in% sequences), ]

# Binarize the dataframe for columns starting with 'r', 't', or 's'
binary_df <- data %>% mutate(across(matches("^(r|t|s)"), ~ ifelse(is.na(.), 0, 1)))

# Sum the binary values per sample
gene_count_per_sample <- colSums(binary_df %>% select(matches("^(r|t|s)")))

# Convert to a dataframe for plotting
gene_count_df <- data.frame(
  sample = names(gene_count_per_sample),
  gene_count = gene_count_per_sample
)

# Create a group column based on the starting letter of the sample
gene_count_df <- gene_count_df %>%
  mutate(group = case_when(
    startsWith(sample, "r") ~ "regular One-tip",
    startsWith(sample, "t") ~ "sorted in Evotips",
    startsWith(sample, "s") ~ "sorted in eppendorfs"
  ))

#create a new column with the number of cells you used per sample
gene_count_df$cell_number = as.numeric(gsub('.*_(\\d+)k.*', '\\1', gene_count_df$sample)) * 1000

# Custom function to remove the labels for the barplot
rename_labels <- function(x) {
  x = sub("^[^_]*_", " ", x) # replace everything before the first underscore
  x = gsub("_"," ", x) # replace underscores with spaces
  x = sub("eppi", "", x)
  
  
}


# Create the bar plot
ggplot(gene_count_df, aes(x = reorder(sample, cell_number), y = gene_count, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Sample type", y = "Number of Proteins") +
  theme_minimal() +
  scale_fill_manual(values = c("regular One-tip" = "#deebf7", "sorted in Evotips" = "#3182bd", "sorted in eppendorfs" = "#9ecae1")) +
  scale_y_continuous(breaks = seq(0, 5000, 1000), limits = c(0, 6000)) +
  scale_x_discrete(labels = rename_labels) + #remane the labels of the plot according to you customized funtion
  geom_text(aes(label = gene_count), nudge_y = -100 , size= 3.5) #annotate protein number in the bar





