# This scripts uses as input the output "data_results" dataframe from proccessing output proteomics files with dep package
# The aim is to benchmark the results from differne softwares used in the proteomics field
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
diann = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/V23_LK/data_results_V23_DIANN_.tsv', stringsAsFactors = F)
pd = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/V23_LK/data_results_v23_plus_minus_biotin_PD.tsv', stringsAsFactors = F)


diann$plus.biotin_vs_minus.biotin_ratio <- as.numeric(diann$plus.biotin_vs_minus.biotin_ratio)

# Convert df$PD to numeric

pd$plus.biotin.lysed.with.NP.40_vs_minus.biotin.lysed.with.NP.40_ratio <- as.numeric(pd$plus.biotin.lysed.with.NP.40_vs_minus.biotin.lysed.with.NP.40_ratio)


### Creating the dfs for downstream visualization

#It finds the common genes present in both data frames based on the name column
intersect_genes = intersect(diann$name, pd$name)

#Filter the diann and pd dataframe to keep only the rows where the name column matches one of the genes in  in intersect_genes.
#This effectively removes any genes in diann that are not present in pd and vice verca
diann_filter = diann[diann$name %in% intersect_genes,]
pd_filter = pd[pd$name %in% intersect_genes,]

#sort on ascending order based on protein name -> super important for the correlation plot
diann_filter = diann_filter[order(diann_filter$name),]
pd_filter = pd_filter[order(pd_filter$name),]


# Extracting only the proteins that are significant based on the logical vector "significant"
diann_true_values <- diann[diann$significant, ]
pd_true_values <- pd[pd$significant, ]


# Extracting overlapping significant proteins based on the logical vector "significant"
intersect_significant_genes = intersect(diann_true_values$name, pd_true_values$name)


#Filter the dfs to keep only the rows where the name column matches one of the genes in  in intersect_genes.
#This effectively removes any genes in diann that are not present in pd and vice verca
diann_significant_filter = diann_true_values[diann_true_values$name %in% intersect_significant_genes,]
pd_significant_filter = pd_true_values[pd_true_values$name %in% intersect_significant_genes,]

#sort on ascending order based on protein name -> super important for the correlation plot
diann_significant_filter = diann_significant_filter[order(diann_significant_filter$name),]
pd_significant_filter = pd_significant_filter[order(pd_significant_filter$name),]


#This line creates new data frames df and df2 containing three columns: gene, DIANN, and PD for the correlation plots
df = data.frame(gene = pd_filter$name, DIANN = diann_filter$plus.biotin_vs_minus.biotin_ratio, PD = pd_filter$plus.biotin.lysed.with.NP.40_vs_minus.biotin.lysed.with.NP.40_ratio,
              PD_significant =  pd_filter$significant, DIANN_significant = diann_filter$significant)
df2 = data.frame(gene = pd_significant_filter$name, DIANN = diann_significant_filter$plus.biotin_vs_minus.biotin_ratio, PD = pd_significant_filter$plus.biotin.lysed.with.NP.40_vs_minus.biotin.lysed.with.NP.40_ratio)


# Create a new column that will identify if a gene is significant in PD or DIANN or both or none so you can plot it in a scatter plot
# it acconts both for -log2foldchange (theshold 1.5) and adj-pvalue (threshold 0.05)
df_scatter <- df %>%
  mutate(Found   = case_when(
    PD_significant == T & DIANN_significant == T  ~ "Both",
    PD_significant== T ~ "PD",
    DIANN_significant == T ~ "DIANN",
    PD_significant == F & DIANN_significant == F ~ "None"
  ))

####Alternative -exclude the proteins that are not identified as significant in either software

df_scatter <- df %>%
  mutate(Found = case_when(
    PD_significant == T & DIANN_significant == T  ~ "Both",
    PD_significant== T ~ "PD",
    DIANN_significant == T ~ "DIANN",
    PD_significant == F & DIANN_significant == F ~ "None"
  )) %>% 
  filter(Sig != "None")


#This line calculates the correlation coefficients (r) between the diann and pd
r = round(cor(df$DIANN, df$PD), digits = 2)
r2 = round(cor(df2$DIANN, df2$PD), digits = 2)


# Identify non-overlapping significant proteins in each dataframe and put them in a separate df 
non_overlapping_diann <- diann_true_values[!(diann_true_values$name %in% pd_true_values$name), ]
non_overlapping_pd <- pd_true_values[!(pd_true_values$name %in% diann_true_values$name), ]

# Sort on ascending order (which is the default in the sort function) based on the p-adj (4rth column of the dfs) 
# which translates to biologically most significant

non_overlapping_diann <- non_overlapping_diann[order(non_overlapping_diann[,4]),, drop = F ]
non_overlapping_pd <- non_overlapping_pd[order(non_overlapping_pd[,4]),, drop = F]


### Visualization

#scatter plot of the correlation of all mutual proteins in DIANN and PD
ggplot(df, aes(x = DIANN, y = PD)) +
    geom_point(color="black", size=1.5) +
    theme_minimal() +
    geom_smooth(method=lm, color="violet", se =F) +
    annotate(geom="text", x= min(df$DIANN) * 0.9, y=max(df$PD) * 0.9, label=paste0('r = ', r),
                     color="black", size = 4) +
    labs(title = "-log2 fold-change correlation plot of mutual proteins in DIANN 1.8.1 and PD", x = "DIANN", y = "PD") +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 10),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0))  # Adjust margin here
  )

#scatter plot of the correlation of all mutual proteins in DIANN and PD with annotations regarding significance 
ggplot(df_scatter, aes(x = DIANN, y = PD, color= Found) ) +
  geom_point( size=1.5) +
  theme_minimal() +
  geom_smooth(method=lm, color="black", se =F) +
  annotate(geom="text", x= min(df$DIANN) * 0.9, y=max(df$PD) * 0.9, label=paste0('r = ', r),
           color= "black", size = 4) +
  labs(title = "-log2 fold-change correlation plot of mutual proteins in DIANN 1.8.1 and PD", x = "DIANN", y = "PD") +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 10),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0))  # Adjust margin here
  )


# Create a venn diagram to visualize the overlap between the output files

venn = venn.diagram(x = list(diann = diann$name, pd = pd$name),
  category.names = c("DIANN", "PD"),
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  main = "Number of proteins that overlap between DIANN 1.8.1 and PD",
  main.cex = 1,  # Adjust main title size
  main.fontfamily = "serif",  # Set font family for the main title
  main.fontface = "bold",  # Set font face for the main title
  main.col = "black",  # Set color for the main title
  sub = "LS174Tw4 cells: Plus-Minus Biotin Strep pull-down (Astral DIA, 30SPD Ionoptics)",
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
ggplot(df2, aes(x = DIANN, y = PD)) +
  geom_point(color="black", size=1.5) +
  theme_minimal() +
  geom_smooth(method=lm, color="violet", se =F) +
  annotate(geom="text", x= min(df2$DIANN) * 0.9, y=max(df2$PD) * 0.9, label=paste0('r = ', r2),
           color="black", size = 4) +
  labs(title = "-log2 fold-change correlation plot of significant proteins in DIANN 1.8.1 and PD", x = "DIANN", y = "PD") +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 10),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0))  # Adjust margin here
  )


# Create a Venn diagram for the overlap of significant values
venn_significant <- venn.diagram(x = list(diann_true_values = diann_true_values$name, pd_true_values = pd_true_values$name),
  category.names = c("DIANN", "PD"),
  filename = NULL,
  output = TRUE ,
  main = "Number of significant proteins that overlap between DIANN 1.8.1 and PD",
  main.cex = 1,  # Adjust main title size
  sub = "LS174T w4 cells: Plus-Minus Biotin Strep pull-down (Astral DIA, 30SPD Ionoptics)", #add subtitle 
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



# Plot a heatmap of the overlapping significant proteins found in DIANN and PD based on p-adj value 
# Assign conditions
condition_1 = "PD" #conditions for heatmaps #Now its like that because we compare the softwares
condition_2 = "DIANN" #conditions for heatmaps # In most cases it would be wt vs mutant etc

#create a new column named "condition" in each df that will contain the conditions for the heatmap
pd_significant_filter <- pd_significant_filter %>%
  mutate(condition = condition_1)

diann_significant_filter <- diann_significant_filter %>%
  mutate(condition = condition_2)


# Extract common column names between the dfs to be merged 
common_columns <- intersect(names(diann_significant_filter), names(pd_significant_filter))

#Merge the 2 dataframes
merged_df <- merge(diann_significant_filter, pd_significant_filter, by = common_columns, all = T)

merged_df = merged_df[order(merged_df$minus.biotin_vs_plus.biotin_p.adj),] 

# Keep only the first 20 rows e.g. the 20 most significant genes
  df_top25 <- merged_df[1:25, ]
  
# Plot a heatmap using ggplot where condition is on the x axis and genes on the y-axis 
# and adjusted p-value are the values
ggplot(df_top25, aes(condition, name, fill = minus.biotin_vs_plus.biotin_p.adj)) + geom_tile() +
labs(
  title = "Top 25 significant proteins identified in both DIANN 1.8.1 and PD",
  subtitle = "Astral DIA (30SPD) Ionoptics",
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


