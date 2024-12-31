library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)

#requires iBAQ info file & data_results file from dep

# Keep only the relevant columns with the iBAQ values
iBAQ_info <- iBAQ_info %>%
  select(Protein.Group, Genes, contains("_iBAQ"))


# Define column positions for assigning the condition -> change each time
cols <- c(3:10)

# Assign conditions based on column positions
df <- iBAQ_info %>%
  pivot_longer(cols = starts_with("F..LK"), names_to = "Sample", values_to = "iBAQ") %>%
  mutate(
    Condition = case_when(
      Sample %in% colnames(iBAQ_info)[cols] ~ "V85_T_B", #condition name -> change each time
      TRUE ~ "other"
    )
  )

# calculate mean iBAQ value for each condition
result <- df %>%
group_by(Genes, Condition) %>%
  summarize(Mean_iBAQ = mean(iBAQ, na.rm = TRUE), .groups = "drop")

# Unpivot the dataframe: Spread Condition into separate columns
unpivoted_result <- result %>%
  pivot_wider(names_from = Condition, values_from = Mean_iBAQ)

#calculate log2 of iBAQ values
unpivoted_result <- unpivoted_result %>%
  mutate(across(where(is.numeric), ~ log2(replace(., . == 0, NA))))


# Merge the dataframes
merged_df <- unpivoted_result %>%
  inner_join(data_results, by = c("Genes" = "name"))


#make groups based on significance & up or downregulation and assign them to a new column

merged_df$Legend <- "Background"
merged_df$Legend[merged_df$V85_T_vs_V85_B_p.adj < 0.05 & merged_df$V85_T_vs_V85_B_ratio > log2(1.5)] <- "Enriched" #for the proteins upregulated on the right -> you can rename it well
merged_df$Legend[merged_df$V85_T_vs_V85_B_p.adj < 0.05 & merged_df$V85_T_vs_V85_B_ratio < -log2(1.5)] <- "Background" #for the proteins downregulated on the left -> you can rename it well

# List of proteins to annotate as hits
hit_genes <- c("SMAD4", "HAND1","TCF3", "TCF4","TCF12","ASCL2", "SMAD3", "SMAD2", "SKI", "SKIL"
               
)  # Add your list of genes here

# Add Hits annotation to the Legend column for genes in hits_genes
merged_df <- merged_df %>%
  mutate(Legend = ifelse(Genes %in% hit_genes, "Hits", Legend))


# Create the scatter plot -> change columns each time
ggplot(merged_df, aes(x = V85_T_B , y = V85_T_vs_V85_B_ratio, color = Legend)) + 
  geom_point() +  # Plot the points
  scale_color_manual(values = c("Enriched" = "#7b3294","Downregulated"="#5ab4ac", "Background" = "azure3", "Hits"= "#ca0020")) + #manual color setting based on Legend column
  theme_minimal() +  # Use a minimal theme
  labs(x = "log2(mean_iBAQ_JD39_neg)", y = "log2(JD39/negative_ctrl)") +# Labels and title
  # Text Annotate the hits
  geom_text_repel(data = subset(merged_df, Genes %in% hit_genes), # Text labels for hits genes
                  aes(label = Genes, fontface = ifelse(Genes %in% hit_genes, "bold", "plain")), 
                  color = "#ca0020", # Color for hits genes
                  nudge_x = 0.1, nudge_y = 0.1, 
                  size = 3) 
  
# Save the merged dataframe as a tsv  file in the specified directory ->always change the file name 
data_folder <- "R:/Group Vermeulen/Lila/Mass_spec_results/V111_LK_HAND1_PD"

write.table(merged_df, file = paste0(data_folder, "/MA_plot_V111.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

