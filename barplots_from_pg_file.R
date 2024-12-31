# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(DT)
library(tibble)
library(shiny)
library(shinydashboard)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(R.utils)
library(forcats)


# Convert columns that match the regex (basically the columns that contain maxLFQ intensities) to numeric
columns_to_convert <- grep("V85|V88|V93", names(data), value = TRUE)
data[columns_to_convert] <- lapply(data[columns_to_convert], as.numeric)


# Function to clean up Gene Names from column "Genes" in DIANN 
data <- data %>%
  mutate(Genes = map(strsplit(Genes, ";"), 1))

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

# Binarize the dataframe for columns starting with what you want -> update each time
binary_df <- data %>% mutate(across(matches("V85|V88|V93"), ~ ifelse(is.na(.), 0, 1)))


# Sum the binary values per sample
#gene_count_per_sample <- colSums(binary_df %>% select(matches("^(|s|r|T)")))

# Sum the binary values per sample
binary_columns <- grep("V85|V88|V93", names(binary_df), value = TRUE)
gene_count_per_sample <- colSums(binary_df[, binary_columns])

# Convert to a dataframe for plotting
gene_count_df <- data.frame(
  sample = names(gene_count_per_sample),
  gene_count = gene_count_per_sample
)

#create an extra column grouping by conditions
gene_count_df <- gene_count_df %>%
  mutate(Category = case_when(
    str_detect(sample, "V85_T") ~ "experiment_1_T_cells",
    str_detect(sample, "V88_T") ~ "experiment_2_T_cells",
    str_detect(sample, "V93_T") ~ "experiment_3_T_cells",
    str_detect(sample, "V85_B") ~ "experiment_1_B_cells",
    str_detect(sample, "V88_B") ~ "experiment_2_B_cells",
    str_detect(sample, "V93_B") ~ "experiment_3_B_cells"
    ))

#IF APPLICABLE: create a new column with the number of cells you used per sample
#gene_count_df$cell_number = as.numeric(gsub('.*_(\\d+)k.*', '\\1', gene_count_df$sample)) * 1000

# Custom function to remove the labels for the barplot
rename_labels <- function(x) {
 # x = sub("^[^_]*_", "", x) # replace everything before the first underscore
  x = gsub("_"," ", x) # replace underscores with spaces
  x = sub("tip_","", x) # replace your desired phrase with nothing
  x = sub(".raw", "", x)
  x = sub("V85","exp_1", x)
  x = sub("V88", "exp_2", x)
  x = sub("V93", "exp_3", x)
  
  
}


#Optional: Rank samples based on the gene_count from highest to lowest
#gene_count_df$sample <- fct_reorder(gene_count_df$sample, gene_count_df$gene_count)


# Create the bar plot
ggplot(gene_count_df, aes(x = sample , y = gene_count, fill = Category)) + #fill based on the groups you created
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "Sample type", y = "Number of Proteins", title = "Sorted in tip") +
  theme_minimal() +
  scale_fill_manual(values = c("experiment_3_B_cells" = "#ffffcc", "experiment_2_B_cells"="#c7e9b4", "experiment_1_B_cells"="#66c2a4", "experiment_3_T_cells"="#edf8fb", "experiment_2_T_cells"="#b3cde3","experiment_1_T_cells"="#8c96c6")) +
  scale_y_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 6000)) + 
  scale_x_discrete(labels = rename_labels) + #remane the labels of the plot according to you customized funtion 
  geom_text(aes(label = gene_count), nudge_y = -120 , size= 3.5) #annotate protein number in the bar 



#Alternative 2: Create the bar plot with reorder in the x axis based on cell number
ggplot(gene_count_df, aes(x = reorder(sample, cell_number), y = gene_count, fill = Category)) + #fill based on the groups you created
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "Sample type", y = "Number of Proteins", title = "Sorted T cells") +
  theme_minimal() +
  scale_fill_manual(values = c("untreated_1k" = "#ffffcc", "untreated_2k"="#c7e9b4", "untreated_5k"="#7fcdbb" )) +
  scale_y_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 5000)) + 
  scale_x_discrete(labels = rename_labels) + #remane the labels of the plot according to you customized funtion
  geom_text(aes(label = gene_count), nudge_y = -128 , size= 3.5) #annotate protein number in the bar



#Alternative 3:
##### If you want to rank the samples explicitly within each Category
# Reorder the samples explicitly within each Category
gene_count_df <- gene_count_df %>%
  group_by(Category) %>%
  mutate(sample_ordered = factor(sample, levels = sample[order(gene_count)])) %>% #gene_count based on descending order or -gene_count based on ascending order
  ungroup()

# Create the bar plot with the reordered sample_ordered variable
ggplot(gene_count_df, aes(x = sample_ordered, y = gene_count, fill = Category)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "Sample type", y = "Number of Proteins", title = "5K cells sorted in tips") +
  theme_minimal() +
  scale_fill_manual(values = c("experiment_3_B_cells" = "#c7e9b4", 
                               "experiment_2_B_cells" = "#c2e699", 
                               "experiment_1_B_cells" = "#66c2a4", 
                               "experiment_3_T_cells" = "#fed98e", 
                               "experiment_2_T_cells" = "#fe9929",
                               "experiment_1_T_cells" = "#d95f0e")) +
  scale_y_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 6000)) + 
  scale_x_discrete(labels = rename_labels) + 
  geom_text(aes(label = gene_count), nudge_y = -120 , size= 3.5)


