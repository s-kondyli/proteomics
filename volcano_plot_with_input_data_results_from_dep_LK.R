library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(purrr)

# Creates a volcano plot using as input data_results table extracted from processing the analyses files with dep package
# In this package the reported adj p-value is plainly reported so you need to calculate the -log10
# The fold-change is already reported as log2fold change
# Clean environment and set seed ----
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)
set.seed(123456)

#Import data from "Import Dataset" and give the object name "data"

#Filter out contaminants
data <- subset(data, !grepl("^KRT*", name))
# Transform column to numeric
data$plus.biotin_vs_minus.biotin_ratio <- as.numeric(data$plus.biotin_vs_minus.biotin_ratio)
data$plus.biotin_vs_minus.biotin_p.adj <- as.numeric(data$plus.biotin_vs_minus.biotin_p.adj)

#select columns for volcano plot, update every time x:ratio, y: p-adj value, third variable the column that contains the genes
df_volcano <- data %>%
  dplyr::select(plus.biotin_vs_minus.biotin_ratio, plus.biotin_vs_minus.biotin_p.adj, name)

colnames(df_volcano) <- c('x','y','Genes')

#make groups based on significance and assign them in a new column
df_volcano$Legend <- "Background"
df_volcano$Legend[df_volcano$y < 0.05 & df_volcano$x > log2(1.5)] <- "Enriched" #here the data on the y-axis is the adj p-value 
df_volcano$Legend[df_volcano$y < 0.05 & df_volcano$x < -log2(1.5)] <- "Enriched" #and on the x-axis the log2(fold-change)
#df_volcano$Legend[data$significant ==T] <- "Significant" #alternative 

# List of proteins to annotate as hits
hit_genes <- c("STK1a1", "FZD1", "FZD2","FZD3","FZD6", "FZD7")  # Add your list of genes here

# Add Hits annotation to the Legend column for genes in hits_genes
df_volcano <- df_volcano %>%
  mutate(Legend = ifelse(Genes %in% hit_genes, "Hits", Legend))

# Visualize
ggplot(df_volcano, aes(x = x, y = -log10(y), color= Legend)) + #color based on the Legend column
  geom_point() +
  labs(x = bquote(Log[2] ( "Fold change")), y = bquote(-Log[10]("adjusted p-value")),
       title = "Differential expression analysis of Colo320",
       subtitle = "BH FDR correction, analyzed with PD",
       caption = "Astral DIA (30SPD)-Ionoptics") +
  scale_color_manual(values = c("Enriched" = "slateblue4", "Background" = "azure3", "Hits"= "red")) + #manual color setting
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), linetype = "dashed") + #lines marking the significance thresholds
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # Text labels only for significant genes
  geom_text_repel(data = subset(df_volcano, Legend == "Enriched"), 
                  aes(label = Genes),
                  nudge_x = 0.1, nudge_y = 0.1, # Adjust the position of labels
                  size = 3) +  # Adjust the size of labels  
  #annotate the hits
  geom_text_repel(data = subset(df_volcano, Genes %in% hit_genes), # Text labels for hits genes
                  aes(label = Genes, fontface = ifelse(Genes %in% hit_genes, "bold", "plain")), 
                  color = "red", # Color for hits genes
                  nudge_x = 0.1, nudge_y = 0.1, 
                  size = 3) + 
  
  # Add text annotations for the conditions -> adjust each time
  annotate("text", x = min(df_volcano$x), y = min(df_volcano$y), label = "untreated", hjust = 0, vjust = 0) + #left
  annotate("text", x = max(df_volcano$x), y = min(df_volcano$y), label = "plus biotin", hjust = 1, vjust = 0) +  # right
  theme_minimal()
