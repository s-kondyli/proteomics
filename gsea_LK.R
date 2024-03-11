library(msigdb)
library(fgsea)
library(tools)
library(ggplot2)
# This script performs GSEA analysis (CC for cellular compartment and BP for bioological processes 
# It uses as input the output table of dep analysis package
# Load data- adjust path each time
data= read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/V18_plus_minus_biotin/data_results_v18_plus_minus_biotin_PD.tsv',
                   stringsAsFactors = F)


msigdb.hs = getMsigdb(org = 'hs', id = 'SYM')
msigdb.hs = appendKEGG(msigdb.hs)

genesets = subsetCollection(msigdb.hs, collection = 'c5' , subcollection = 'GO:CC') 
genesets = genesets[grepl('^GOCC', names(genesets))] # CC for cellular compartment and BP for biological processes
pathways = geneIds(genesets)

# rank genes based on the -log2fold ratio -> adjust each time
ranked_genes = data$plus.biotin_vs_minus.biotin_ratio
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


df = data.frame(pathway = fgseaRes$pathway, padj = fgseaRes$padj, NES = fgseaRes$NES, cat = 'CC')


# New df -> Combine the top x rows with the highest NES values 
# and the top x rows with the lowest NES values from df -> you can adjust it

df_plot = rbind(df[order(df$NES, decreasing = T),][1:20,], df[order(df$NES),][1:5,]) 


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
  title = "GSEA analysis for surface biotinylation",
  subtitle = "Colo 320 cells +/- biotin lysed with the commercial lysis buffer",
  caption = "Astral DIA (30SPD) Ionoptics analyzed with PD", )



