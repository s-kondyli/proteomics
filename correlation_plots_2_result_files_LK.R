# This scripts uses as input the output "data_results" dataframe from proccessing output proteomics files with dep package
# The aim is to benchmark the results from different softwares used in the proteomics field
# Possibilities: It creates correlation plots of the log2foldchange of the common proteins of the softwares
# venn diagramms and heatmaps


library(ggplot2)
library(reshape2)
library(VennDiagram)
library(gridExtra)
library(RColorBrewer)
library(dplyr)
library(tidyr)

### Load data

#read the dataframes from the saved folder -> adjust each time
result_df_1 = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/V88_PBMCs/data_results_fdrtool_V88_PBMCs_comma_DIANN.tsv', stringsAsFactors = F)
result_df_2 = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/V88_PBMCs/data_results_fdrtool_V88_PBMCs_period_DIANN.tsv', stringsAsFactors = F)


result_df_1$T_cells_vs_B_cells._ratio <- as.numeric(result_df_1$T_cells_vs_B_cells._ratio)

# Convert df$result_df_2 to numeric

result_df_2$T_cells_vs_B_cells._ratio <- as.numeric(result_df_2$T_cells_vs_B_cells._ratio)


### Creating the dfs for downstream visualization

#It finds the common genes present in both data frames based on the name column
intersect_genes = intersect(result_df_1$name, result_df_2$name)

#Filter the result_df_1 and result_df_2 dataframe to keep only the rows where the name column matches one of the genes in intersect_genes.
#This effectively removes any genes in result_df_1 that are not present in result_df_2 and vice verca
result_df_1_filter = result_df_1[result_df_1$name %in% intersect_genes,]
result_df_2_filter = result_df_2[result_df_2$name %in% intersect_genes,]

#sort on ascending order based on protein name -> super important for the correlation plot
result_df_1_filter = result_df_1_filter[order(result_df_1_filter$name),]
result_df_2_filter = result_df_2_filter[order(result_df_2_filter$name),]


# Extracting only the proteins that are significant based on the logical vector "significant"
result_df_1_true_values <- result_df_1[result_df_1$significant, ]
result_df_2_true_values <- result_df_2[result_df_2$significant, ]


# Extracting overlapping significant proteins based on the logical vector "significant"
intersect_significant_genes = intersect(result_df_1_true_values$name, result_df_2_true_values$name)


#Filter the dfs to keep only the rows where the name column matches one of the genes in  in intersect_genes.
#This effectively removes any genes in result_df_1 that are not present in result_df_2 and vice verca
result_df_1_significant_filter = result_df_1_true_values[result_df_1_true_values$name %in% intersect_significant_genes,]
result_df_2_significant_filter = result_df_2_true_values[result_df_2_true_values$name %in% intersect_significant_genes,]

#sort on ascending order based on protein name -> super important for the correlation plot
result_df_1_significant_filter = result_df_1_significant_filter[order(result_df_1_significant_filter$name),]
result_df_2_significant_filter = result_df_2_significant_filter[order(result_df_2_significant_filter$name),]


#This line creates new data frames df and df2 containing three columns: gene, result_df_1, and result_df_2 for the correlation plots
df = data.frame(gene = result_df_2_filter$name, result_df_1 = result_df_1_filter$T_cells_vs_B_cells._ratio, result_df_2 = result_df_2_filter$T_cells_vs_B_cells._ratio,
                result_df_2_significant =  result_df_2_filter$significant, result_df_1_significant = result_df_1_filter$significant)
df2 = data.frame(gene = result_df_2_significant_filter$name, result_df_1 = result_df_1_significant_filter$T_cells_vs_B_cells._ratio, result_df_2 = result_df_2_significant_filter$T_cells_vs_B_cells._ratio)


# Create a new column that will identify if a gene is significant in result_df_2 or result_df_1 or both or none so you can plot it in a scatter plot
# it acconts both for -log2foldchange (theshold 1.5) and adj-pvalue (threshold 0.05)
df_scatter <- df %>%
  mutate(Found   = case_when(
    result_df_2_significant == T & result_df_1_significant == T  ~ "Both",
    result_df_2_significant== T ~ "result_df_2",
    result_df_1_significant == T ~ "result_df_1",
    result_df_2_significant == F & result_df_1_significant == F ~ "None"
  ))

####Alternative -exclude the proteins that are not identified as significant in either software

df_scatter <- df %>%
  mutate(Found = case_when(
    result_df_2_significant == T & result_df_1_significant == T  ~ "Both",
    result_df_2_significant== T ~ "result_df_2",
    result_df_1_significant == T ~ "result_df_1",
    result_df_2_significant == F & result_df_1_significant == F ~ "None"
  )) %>% 
  filter(Sig != "None")


#This line calculates the correlation coefficients (r) between the result_df_1 and result_df_2
r = round(cor(df$result_df_1, df$result_df_2), digits = 2)
r2 = round(cor(df2$result_df_1, df2$result_df_2), digits = 2)


# Identify non-overlapping significant proteins in each dataframe and put them in a separate df 
non_overlapping_result_df_1 <- result_df_1_true_values[!(result_df_1_true_values$name %in% result_df_2_true_values$name), ]
non_overlapping_result_df_2 <- result_df_2_true_values[!(result_df_2_true_values$name %in% result_df_1_true_values$name), ]

# Sort on ascending order (which is the default in the sort function) based on the p-adj (4rth column of the dfs) 
# which translates to biologically most significant

non_overlapping_result_df_1 <- non_overlapping_result_df_1[order(non_overlapping_result_df_1[,4]),, drop = F ]
non_overlapping_result_df_2 <- non_overlapping_result_df_2[order(non_overlapping_result_df_2[,4]),, drop = F]


### Visualization

#scatter plot of the correlation of all mutual proteins in result_df_1 and result_df_2
ggplot(df, aes(x = result_df_1, y = result_df_2)) +
  geom_point(color="black", size=1.5) +
  theme_minimal() +
  geom_smooth(method=lm, color="violet", se =F) +
  annotate(geom="text", x= min(df$result_df_1) * 0.9, y=max(df$result_df_2) * 0.9, label=paste0('r = ', r),
           color="black", size = 4) +
  labs(title = "-log2 fold-change correlation plot of mutual proteins in result_df_1 1.8.1 and result_df_2", x = "result_df_1", y = "result_df_2") +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 10),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0))  # Adjust margin here
  )

#scatter plot of the correlation of all mutual proteins in result_df_1 and result_df_2 with annotations regarding significance 
ggplot(df_scatter, aes(x = result_df_1, y = result_df_2, color= Found) ) +
  geom_point( size=1.5) +
  theme_minimal() +
  geom_smooth(method=lm, color="black", se =F) +
  annotate(geom="text", x= min(df$result_df_1) * 0.9, y=max(df$result_df_2) * 0.9, label=paste0('r = ', r),
           color= "black", size = 4) +
  labs(title = "-log2 fold-change correlation plot of mutual proteins in result_df_1 1.8.1 and result_df_2", x = "result_df_1", y = "result_df_2") +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 10),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0))  # Adjust margin here
  )


# Create a venn diagram to visualize the overlap between the output files

venn = venn.diagram(x = list(result_df_1 = result_df_1$name, result_df_2 = result_df_2$name),
                    category.names = c("result_df_1", "result_df_2"),
                    filename = NULL,
                    output = TRUE,
                    imagetype = "png",
                    main = "Number of proteins that overlap between result_df_1 1.8.1 and result_df_2",
                    main.cex = 1,  # Adjust main title size
                    main.fontfamily = "serif",  # Set font family for the main title
                    main.fontface = "bold",  # Set font face for the main title
                    main.col = "black",  # Set color for the main title
                    sub = "LS174Tw4 cells: Plus-Minus Biotin Strep pull-down (Astral DIA, 30Sresult_df_2 Ionoptics)",
                    sub.col = "blueviolet",  # Set color for the subtitle
                    sub.cex = 0.7,  # Adjust subtitle size
                    sub.fontfamily = "sans",  # Set font family for the subtitle
                    sub.fontface = "italic",  # Set font face for the subtitle
                    category.cex = 1.5,  # Adjust category names size
                    category.fontfamily = "serif",  # Set font family for category names
                    category.fontface = "bold",  # Set font face for category names
                    category.col = c("blue", "lightgreen"))  # Set colors for the categories

# Plot the Venn diagram
grid.draw(venn)

### Significant values
#scatter plot of the correlation for significant proteins
ggplot(df2, aes(x = result_df_1, y = result_df_2)) +
  geom_point(color="black", size=1.5) +
  theme_minimal() +
  geom_smooth(method=lm, color="violet", se =F) +
  annotate(geom="text", x= min(df2$result_df_1) * 0.9, y=max(df2$result_df_2) * 0.9, label=paste0('r = ', r2),
           color="black", size = 4) +
  labs(title = "-log2 fold-change correlation plot of significant proteins in result_df_1 1.8.1 and result_df_2", x = "result_df_1", y = "result_df_2") +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 10),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0))  # Adjust margin here
  )


# Create a Venn diagram for the overlap of significant values
venn_significant <- venn.diagram(x = list(result_df_1_true_values = result_df_1_true_values$name, result_df_2_true_values = result_df_2_true_values$name),
                                 category.names = c("result_df_1", "result_df_2"),
                                 filename = NULL,
                                 output = TRUE ,
                                 main = "Number of significant proteins that overlap between result_df_1 1.8.1 and result_df_2",
                                 main.cex = 1,  # Adjust main title size
                                 sub = "LS174T w4 cells: Plus-Minus Biotin Strep pull-down (Astral DIA, 30Sresult_df_2 Ionoptics)", #add subtitle 
                                 main.fontfamily = "serif",  # Set font family for the main title
                                 main.fontface = "bold",  # Set font face for the main title
                                 main.col = "black",  # Set color for the main title
                                 sub.col = "blueviolet",  # Set color for the subtitle
                                 sub.cex = 0.7,  # Adjust subtitle size
                                 sub.fontfamily = "sans",  # Set font family for the subtitle
                                 sub.fontface = "italic",  # Set font face for the subtitle
                                 category.cex = 1.5,  # Adjust category names size
                                 category.fontfamily = "serif",  # Set font family for category names
                                 category.fontface = "bold",  # Set font face for category names
                                 category.col = c("blue", "lightgreen") ) # Set colors for the categories

grid.draw(venn_significant)



# Plot a heatmap of the overlapping significant proteins found in result_df_1 and result_df_2 based on p-adj value 
# Assign conditions
condition_1 = "result_df_2" #conditions for heatmaps #Now its like that because we compare the softwares
condition_2 = "result_df_1" #conditions for heatmaps # In most cases it would be wt vs mutant etc

#create a new column named "condition" in each df that will contain the conditions for the heatmap
result_df_2_significant_filter <- result_df_2_significant_filter %>%
  mutate(condition = condition_1)

result_df_1_significant_filter <- result_df_1_significant_filter %>%
  mutate(condition = condition_2)


# Extract common column names between the dfs to be merged 
common_columns <- intersect(names(result_df_1_significant_filter), names(result_df_2_significant_filter))

#Merge the 2 dataframes
merged_df <- merge(result_df_1_significant_filter, result_df_2_significant_filter, by = common_columns, all = T)

merged_df = merged_df[order(merged_df$minus.biotin_vs_plus.biotin_p.adj),] 

# Keep only the first 20 rows e.g. the 20 most significant genes
df_top25 <- merged_df[1:25, ]

# Plot a heatmap using ggplot where condition is on the x axis and genes on the y-axis 
# and adjusted p-value are the values
ggplot(df_top25, aes(condition, name, fill = minus.biotin_vs_plus.biotin_p.adj)) + geom_tile() +
  labs(
    title = "Top 25 significant proteins identified in both result_df_1 1.8.1 and result_df_2",
    subtitle = "Astral DIA (30Sresult_df_2) Ionoptics",
    x = "Conditions",
    y = "Enriched Proteins",
    caption = "Source: LS174Tw4 cells +/- biotin treatment followed by strep pull-down",
    fill = "adjusted p-value" #legend name
  ) +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 11),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 2, b = 1, l = 1))  # Adjust margin here
  )


