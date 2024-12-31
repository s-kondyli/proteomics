library(msigdb)
library(fgsea)
library(tools)
library(ggplot2)
library(ReactomePA)

# Load data- adjust path each time
data= read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all/all_filtered_2_unique_peptides/data_PBMCs_all_BH_correction_DIANN_1.9.2.tsv',
                 stringsAsFactors = F)


msigdb.hs = getMsigdb(org = 'hs', id = 'SYM')
msigdb.hs = appendKEGG(msigdb.hs)


genesets = subsetCollection(msigdb.hs, collection = 'c7' , subcollection = 'GO:BP') 
genesets = genesets[grepl('^GOBP', names(genesets))] ## CC for cellular compartment and BP for biological processes, MF for molecular functions
pathways = geneIds(genesets)

# rank genes based on the -log2fold ratio -> adjust each time
ranked_genes = data$V93_T_vs_V93_B_ratio
names(ranked_genes) = data$name
ranked_genes = sort(ranked_genes, decreasing = T)


# creates a df GSEA results. Each row corresponds to a tested pathway
fgseaRes = fgsea(pathways = pathways,
                 stats    = ranked_genes,
                 minSize  = 15,
                 maxSize  = 500,
                 eps = 0)


fgseaRes = fgseaRes[order(fgseaRes$padj),] #order on ascending order based on the p-adj so biologically more significant

# Create dataframe with enriched pathways, adjusted p-values and NES:Normalized Enrichment Scores


df = data.frame(pathway = fgseaRes$pathway, padj = fgseaRes$padj, NES = fgseaRes$NES, cat = 'BP') #change the cat to match the gsea you do e.g. BP, CC or MF


# New df -> Combine the top 10 rows with the highest NES values 
# and the top 10 rows with the lowest NES values from df -> you can adjust it

df_plot = rbind(df[order(df$NES, decreasing = T),][1:10,], df[order(df$NES),][1:10,]) 


#Remove prefixes matching the pattern 'GO.{2}_'.
#Replace underscores with spaces.
#Convert the text to lowercase.
#Convert the text to title case using a custom function (toTitleCase).

df_plot$pathway = sapply(df_plot$pathway, function(x){
  
  p = gsub('GO.{2}_', '', x) 
  p = gsub('_', ' ', p)
  p = tolower(p)
  p = toTitleCase(p)
  return(p)
})

### Plot it 
ggplot(df_plot, aes(x = cat, y = reorder(pathway, NES), size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = '#3300CC', mid = 'white', high = '#FF0000') +
  theme_minimal() +
  theme(axis.title = element_blank()) +
  labs(
    title = "GSEA analysis",
    subtitle = "T vs B cells from PBMCs-3",
    caption = "Astral DIA (30SPD) Ionoptics analyzed with DIANN 1.9.2", )

###### Venn diagram to see overlap of genes between different pathways ONLY FOR SIGNIFICANT PATHWAYS

# Drop the rows that have a padj >= 0.05, so only keep significant ones

fgseaRes_significant_only <- subset(fgseaRes, padj < 0.05)

# Find pathways of interest -> adjust each time
mitochondrial_matching_pathways <- fgseaRes_significant_only$pathway[grep("^GOCC_MITOCHONDRIAL_.", fgseaRes_significant_only$pathway)]
cell_surface_matching_pathways <- fgseaRes_significant_only$pathway[grep("^GOCC.*(?:PLASMA|SURFACE)", fgseaRes_significant_only$pathway)]
ribosomal_matching_pathways <- fgseaRes_significant_only$pathway[grep("^GOCC.*(?:RIBOSO)", fgseaRes_significant_only$pathway)]

# Extract genes associated with matching pathways
mitochondrial_matching_genes <- unlist(fgseaRes_significant_only$leadingEdge[match(mitochondrial_matching_pathways, fgseaRes_significant_only$pathway)])
cell_surface_matching_genes <- unlist(fgseaRes_significant_only$leadingEdge[match(cell_surface_matching_pathways, fgseaRes_significant_only$pathway)])
ribosomal_matching_genes <- unlist(fgseaRes_significant_only$leadingEdge[match(ribosomal_matching_pathways, fgseaRes_significant_only$pathway)])

# plots the venn diagram -> adjust the sets each time
venn.plot <- venn.diagram(
  x = list(set1 = ribosomal_matching_genes, set2 = cell_surface_matching_genes),
  category.names = c("ribosomal", "mitochondrial"),
  filename = NULL
)
grid.draw(venn.plot)



