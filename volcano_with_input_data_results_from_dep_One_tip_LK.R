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

#Import data from "Import Dataset"

# Transform column to numeric
data$V85_T_vs_V85_B_p.adj<- as.numeric(data$V85_T_vs_V85_B_p.adj)
data$V85_T_vs_V85_B_ratio <- as.numeric(data$V85_T_vs_V85_B_ratio)

#select columns for volcano plot, update every time x:ratio, y: p-adj value
df_volcano <- data %>%
  dplyr::select(V85_T_vs_V85_B_ratio, V85_T_vs_V85_B_p.adj, name)

colnames(df_volcano) <- c('x','y','Genes')

#make groups based on significance & up or downregulation and assign them in a new column
#here the data on the y-axis is the adj p-value 
#and on the x-axis the log2(fold-change)
df_volcano$Legend <- "Background"
df_volcano$Legend[df_volcano$y < 0.05 & df_volcano$x > log2(1.5)] <- "T_cells" #for the proteins upregulated on the right -> you can rename it well
df_volcano$Legend[df_volcano$y < 0.05 & df_volcano$x < -log2(1.5)] <- "B_cells" #for the proteins downregulated on the left -> you can rename it well
#df_volcano$Legend[data$significant ==T] <- "Significant" #alternative 

# List of proteins to annotate as hits
hit_genes <- c("CD3E", "CD3D","CD247", "CD8A","CD28", "CD5","IL7R","PTPRC","SELL", "CD4" 
               , "CD19", "MS4A1", "CR2", "CD22" , "TCL1A", "IGHM","CD74" , "CD79A", "CD79B", "CD3G", "CD8B", "TCF7"
               )  # Add your list of genes here

# Add Hits annotation to the Legend column for genes in hits_genes
df_volcano <- df_volcano %>%
  mutate(Legend = ifelse(Genes %in% hit_genes, "Hits", Legend))



ggplot(df_volcano, aes(x = x, y = -log10(y), color= Legend)) + #color based on the Legend column
  geom_point() +
  labs(x = bquote(Log[2] ( "Fold change")), y = bquote(-Log[10]("adjusted p-value")),
       title = "Differential expression analysis of T vs B cells sorted from PBMCs-1'",
       subtitle = "BH correction, analyzed with DIANN 1.9.2",
       caption = "Astral DIA (30SPD)-Ionoptics") +
  scale_color_manual(values = c("T_cells" = "#7b3294","B_cells"="#5ab4ac", "Background" = "azure3", "Hits"= "red")) + #manual color setting based on Legend column
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), linetype = "dashed") + #lines marking the significance thresholds
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # Text annotate T_cells proteins -> if you don't want text annotations you can skip this part
 # geom_text_repel(data = subset(df_volcano, Legend == "T_cells"), 
 #                 aes(label = Genes),
 #                 nudge_x = 0.1, nudge_y = 0.1, # Adjust the position of labels
  #                size = 3) +  # Adjust the size of labels  
  
  
  # Text Annotate B_cells proteins -> if you don't want text annotations you can skip this part
  #geom_text_repel(data = subset(df_volcano, Legend == "B_cells"), 
  #                aes(label = Genes),
   #               nudge_x = -0.1, nudge_y = 0.1,
   #               size = 3) +
  
  # Text Annotate the hits
  geom_text_repel(data = subset(df_volcano, Genes %in% hit_genes), # Text labels for hits genes
                  aes(label = Genes, fontface = ifelse(Genes %in% hit_genes, "bold", "plain")), 
                  color = "red", # Color for hits genes
                  nudge_x = 0.1, nudge_y = 0.1, 
                  size = 3) + 
  
  # Add text annotations for the conditions -> adjust each time
  annotate("text", x = min(df_volcano$x), y = min(df_volcano$y), label = "B cells", hjust = 0, vjust = 0) + #left
  annotate("text", x = max(df_volcano$x), y = min(df_volcano$y), label = "T cells", hjust = 1, vjust = 0) +  # right
  theme_minimal()
