library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(purrr)

# Creates a volcano plot using as input data_results table extracted from processing the analyses files with dep package
# In this package the reported adj p-value is plainly reported so you need to calculate the -log10
# The fold-change is already reported as log2fold change

#Import dat from "Import Dataset"

# Transform column to numeric
data$wt_vs_PDL2_ko_p.adj <- as.numeric(data$wt_vs_PDL2_ko_p.adj)
data$wt_vs_PDL2_ko_ratio <- as.numeric(data$wt_vs_PDL2_ko_ratio)

#select columns for volcano plot, update every time
df_volcano <- data %>%
  dplyr::select(wt_vs_PDL2_ko_ratio, wt_vs_PDL2_ko_p.adj, name)

colnames(df_volcano) <- c('x','y','Genes')

#make groups based on significance and assign them in a new column
df_volcano$Legend <- "Background"
df_volcano$Legend[df_volcano$y < 0.05 & df_volcano$x > log2(1.5)] <- "Significant" #here the data on the y-axis is the adj p-value 
df_volcano$Legend[df_volcano$y < 0.05 & df_volcano$x < -log2(1.5)] <- "Significant" #and on the x-axis the log2(fold-change)
#df_volcano$Legend[data$significant ==T] <- "Significant" #alternative 

# List of proteins to annotate as hits
hits_genes <- c("PDCD1LG2", "Gene10", "Gene20")  # Add your list of genes here


ggplot(df_volcano, aes(x = x, y = -log10(y), color= Legend)) + #color based on the Legend column
  geom_point() +
  labs(x = bquote(Log[2] ( "Fold change")), y = bquote(-Log[10]("adjusted p-value")),
       title = "Differentially expressed proteins",
       subtitle = "Benjamini-Hochberg FDR control-DIANN 1.8.1 data",
       caption = "Astral DIA (30SPD)-Ionoptics") +
  scale_color_manual(values = c("Significant" = "darkorchid", "Background" = "azure3")) + #manual color setting
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), linetype = "dashed") + #lines marking the significance thresholds
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # Text labels only for significant genes
  geom_text_repel(data = subset(df_volcano, Legend == "Significant"), 
                  aes(label = Genes),
                  nudge_x = 0.1, nudge_y = 0.1, # Adjust the position of labels
                  size = 3) +  # Adjust the size of labels  
  #annotate the hits
  geom_text_repel(data = subset(df_volcano, Genes %in% hits_genes ), # Text labels for hits genes
                  aes(label = Genes, fontface = ifelse(Genes %in% hits_genes, "bold", "plain")), 
                  color = "#333000", # Color for hits genes
                  nudge_x = 0.1, nudge_y = 0.1, 
                  size = 3) + 
  
  # Add text annotations for the conditions -> adjust each time
  annotate("text", x = min(df_volcano$x), y = min(df_volcano$y), label = "PDL2 KO", hjust = 0, vjust = 0) + #left
  annotate("text", x = max(df_volcano$x), y = min(df_volcano$y), label = "wt", hjust = 1, vjust = 0) + #right
  theme_minimal()
