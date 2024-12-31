library(DEP)
library(ggplot2)
library(msigdb)
library(fgsea)
library(tools)
library(SummarizedExperiment)
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)
library(RcisTarget)
library(cowplot)

# Load Data

data("motifAnnotations_hgnc_v9") #Needed for motif analysis
motifAnnotations = motifAnnotations_hgnc_v9

# Prepare Data

prepare_se = function(pg_matrix, expDesign, missing_thr = 0, impute = T){

  pg_uniq = make_unique(pg_matrix, 'Genes', 'Protein.Group', delim = ';')
  pat = paste(expDesign$label, collapse = '|')
  lfq = grep(pat, colnames(pg_uniq))
  se = make_se(pg_uniq, lfq, expDesign)
  se_filt = filter_missval(se, missing_thr)
  se = normalize_vsn(se_filt)
  p = plot_normalization(se_filt, se)
  grid.newpage()
  grid.draw(p)

  if(impute){se = impute(se, 'MinProb')}

  return(se)
}

# Peptide level

summarize_peptide_intensities = function(pr_matrix, id_column = 'Protein.Group', peptide_level = 'Stripped.Sequence'){
  # Helper function for other peptide functions

  ints = as.matrix(pr_matrix[,11:ncol(pr_matrix)])
  ints[is.na(ints)] = 0
  ints = aggregate(ints, list(peptide = pr_matrix[,peptide_level]), 'sum')
  rownames(ints) = ints[,1]
  ints = as.matrix(ints[,-1])
  ints[ints == 0] = NA
  ints = as.data.frame(ints)

  return(ints)
}



get_nPep_prMatrix = function(pr_matrix, id_column = 'Protein.Group', peptide_level = 'Stripped.Sequence'){

  ints = summarize_peptide_intensities(pr_matrix, id_column, peptide_level)

  ints$protein = pr_matrix[match(rownames(ints), pr_matrix[,peptide_level]), id_column]

  npep = aggregate(ints[,-ncol(ints)], list(protein = ints$protein), FUN = function(x){length(na.omit(x))})
  n_total = tapply(rownames(ints), ints$protein, length)
  npep$n_total = n_total[match(npep$protein, names(n_total))]


  #pr_matrixs = merge(pr_matrix, npep, by.x = id_column, by.y = 'protein', suffixes = c('_intensity', '_npep'))

  return(npep)
}



add_peptide_numbers = function(pg_matrix, peptide_numbers, id_column = 'Protein.Group'){

  pg_matrix = merge(pg_matrix, peptide_numbers, by.x = id_column, by.y = 'protein', suffixes = c('_intensity', '_npep'))
  return(pg_matrix)
}



get_median_intensities_prMatrix = function(pr_matrix, id_column = 'Protein.Group', peptide_level = 'Stripped.Sequence'){

  ints = summarize_peptide_intensities(pr_matrix, id_column, peptide_level)
  ints$protein = pr_matrix[match(rownames(ints), pr_matrix[,peptide_level]), id_column]

  pep_median = aggregate(ints[,-ncol(ints)], list(protein = ints$protein), FUN = function(x){median(x, rm.na = T)})
  return(pep_median)
}

plot_peptide_intensities = function(rep, protein, sequence){

  peps = rep[grep(protein, rep$Genes),]
  peps_wide = acast(peps, Stripped.Sequence ~ Run, fun.aggregate = sum, value.var = 'Precursor.Normalised')
  peps_wide = peps_wide[!rowSums(peps_wide) == 0,]

  locations = sapply(rownames(peps_wide), function(x){

    loc = gregexpr(x, sequence)[[1]]
    start = loc[1]
    end = start + (attr(loc, 'match.length') - 1)

    return(c(start = start, end = end))
  })

  peps_wide = cbind(t(locations), log2(peps_wide))
  peps_wide = as.data.frame(peps_wide[order(peps_wide[,1]),])
  peps_wide$sequence = rownames(peps_wide)

  data_long = melt(peps_wide, id.vars = c('sequence', 'start', 'end'), variable.name = 'sample', value.name = 'intensity')
  data_long$condition = gsub('_\\d$', '', data_long$sample)
  data_long$intensity[is.infinite(data_long$intensity)] = 0.01
  data_long$pos = paste0('AA: ', data_long$start, '-', data_long$end)
  data_long$sequence = reorder(data_long$sequence, data_long$start)
  data_long$pos = reorder(data_long$pos, data_long$start)

  ggplot(data_long, aes(x = sample, y = intensity, fill = condition)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    facet_wrap(~ pos) +
    theme_classic() +
    labs(x = '', y = 'log2(Peptide Intensity)') +
    theme(axis.text.x = element_blank(),
          strip.text = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 18))
}



# Analyze DEP

recode_sig_col = function(res, pval_cutoff = 0.05, fc_cutoff = 1){
  # helper function for get_DEPresults

  data = res
  ratio_col = grep('ratio', colnames(data), value = T)
  padj_col = grep('bh', colnames(data), value = T)
  sig_col = grep('sig', colnames(data), value = T)[1]
  data[,sig_col] = ifelse(abs(data[,ratio_col]) > fc_cutoff & data[,padj_col] < pval_cutoff, T, F)
  return(data)
}



get_DEPresults = function(se, condition1 = NULL, condition2 = NULL, alpha = 0.05, lfc = 1, type = 'manual', fdr_type = 'BH'){


  if (type == 'manual'){

    pat = paste(condition1, condition2, sep = '|')
    se = se[,grep(pat, colnames(se))]
    test = paste0(condition1, '_vs_', condition2)
    dep = test_diff(se, type = 'manual', test = test)

    dep <- add_rejections(dep, alpha = alpha, lfc = lfc)
    res = get_results(dep)
    res = res[order(res$name),]

    res$bh = p.adjust(res[,grep('p.val', colnames(res))], method = 'BH')


    fc_col = grep('ratio', colnames(res))
    rd = as.data.frame(rowData(dep))
    diff_col = grep('diff', colnames(rd))
    res[,fc_col] = rd[,diff_col]

    if (fdr_type == 'BH'){
      res[,grep('p.adj', colnames(res))] = res$bh
      res = recode_sig_col(res, alpha, lfc)

    }
  }

  return(res)
}



perform_PCA = function(se, top_n = 500){

  data = assay(se)
  vars = rowVars(t(scale(t(data), scale = F)))
  data = data[order(vars, decreasing = T),][1:top_n,]

  pca = prcomp(t(data))
  pcs = pca$x

  loadings = pca$rotation
  return(list(pcs = pcs, loadings = loadings, summary = summary(pca)))
}

## Plot Functions

add_standardTheme  = function(plot){

  plot + theme(axis.text = element_text(color = 'black', size = 18),
            axis.title = element_text(color = 'black', size = 18))
}



plotVolcano = function(data, pval_cutoff = 0.05, fc_cutoff = 1, label = 'sig', top_n = NULL){

  ratio_col = grep('ratio', colnames(data), value = T)
  padj_col = grep('p.adj', colnames(data), value = T)
  #padj_col = grep('bh', colnames(data), value = T)
  sig_col = grep('sig', colnames(data), value = T)[1]

  data[,sig_col] = ifelse(abs(data[,ratio_col]) > fc_cutoff & data[,padj_col] < pval_cutoff, T, F)


  if (length(label) == 1 & 'sig' %in% label){
    data$label = ifelse(data[[sig_col]], data$name, '')

    if (!is.null(top_n)){
      data_sig = data[data[[sig_col]],]
      data_sig = data_sig[order(data_sig[,padj_col],decreasing = F),]
      sig_pos = data_sig[data_sig[[ratio_col]] > 0,]
      sig_neg = data_sig[data_sig[[ratio_col]] < 0,]

      labels_to_keep = c(sig_pos[1:top_n, 'name'], sig_neg[1:top_n,'name'])
      #labels_to_keep = data_sig[c(1:top_n, (nrow(data_sig) - (top_n -1)):nrow(data_sig)),'name']
      data$label = ifelse(data$name %in% labels_to_keep, data$label, '')
    }
  }
  else{
    data$label = ifelse(data[,'name'] %in% label, data[,'name'], '')
  }

  title = gsub('_p.val', '', colnames(data)[4])
  data$col = ifelse(data$label != '', 'red3', 'grey70')
  #data$col = factor(data$col, levels = c('red3', 'grey70'))
  data = data[order(data$col),]

  p = ggplot(data, aes(x = .data[[ratio_col]], y = -log10(.data[[padj_col]]), label = label)) +
    geom_point(color = data$col) +
    geom_hline(yintercept = -log10(pval_cutoff), color = 'red', linetype = 'dashed', linewidth = 1) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), color = 'red', linetype = 'dashed', linewidth = 1) +
    geom_text_repel(max.overlaps = Inf, min.segment.length = 0.01) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 14))

  p = add_standardTheme(p)
  return(p)

}



get_common_genes = function(res1, res2, ..., pval_cutoff = 0.05, fc_cutoff = 1, type = 'intersect'){

  res_list = list(res1, res2, ...)

  fc_columns  = sapply(res_list, function(x){grep('ratio', colnames(x), value = T)})
  padj_columns = sapply(res_list, function(x){grep('p.adj', colnames(x), value = T)})


  gois = lapply(1:length(res_list), function(x){
    fc = fc_columns[x]
    padj = padj_columns[x]
    res = res_list[[x]]
    genes = res[res[,padj] < pval_cutoff & abs(res[,fc]) > fc_cutoff,'name']
  })

  if (type == 'intersect'){
    return(Reduce(intersect, gois))
  }
  else{
    return(Reduce(union, gois))
  }
}



plot_venn_diagram = function(res1, res2, ..., colors = NULL, names = NULL){

  l = list(res1, res2, ...)
  genes = lapply(l, function(x){
    x[x$significant, 'name']
  })

  if (is.null(names)){names(genes) = paste0('res', 1:length(l))}
  else {names(genes) = names}

  if (is.null(colors)){
    if(length(l) > 2){col_vals = brewer.pal(length(l), 'Set2')}
    else {col_vals = brewer.pal(3, 'Set2')[1:length(l)]}
  }
  else {
    col_vals = colors
  }

  vd = venn.diagram(genes, NULL, T, fill = col_vals)
  grid.newpage()
  grid.draw(vd)
}



plot_nDEP_barplot = function(res1, res2, ..., names = NULL, res_list = NULL){

  if (is.null(res_list)){l = list(res1, res2, ...)}
  else {l = res_list}

  genes = lapply(l, function(x){
    #x = recode_sig_col(x)
    x[x$significant, 'name']
  })

  if (is.null(names)){names(genes) = paste0('res', 1:length(l))}
  else {names(genes) = names}

  df = data.frame(Group = names(genes),
                  nProts = lengths(genes))

  if(!is.null(names)){df$Group = factor(df$Group, levels = names)}

  p = ggplot(df, aes(x = Group, y = nProts, fill = Group)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.75) +
    theme_classic() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 21)) +
    labs(x = '', y = 'No. Significant Proteins') +
    scale_fill_brewer(palette = 'Set2') +
    theme(aspect.ratio = 2 / (nrow(df) - 1))

  p = add_standardTheme(p)
  return(p)
}

## Gene set enrichment analysis

load_msigdb = function(){
  msigdb.hs = getMsigdb(org = 'hs', id = 'SYM')
  msigdb.hs = appendKEGG(msigdb.hs)
  return(msigdb.hs)
}



get_genesets = function(db, collection = c('c2', 'c5', 'h'), subcollection = c('GOBP', 'GOCC', 'GOMF', 'KEGG')){


  genesets = subsetCollection(db, collection = collection)
  if(collection != 'h'){
    pat = paste0('^', subcollection)
    genesets = genesets[grepl(pat, names(genesets))]
  }
  pathways = geneIds(genesets)
}



get_ranked_genes = function(res){

  fc_col = grep('ratio', colnames(res), value = T)
  padj_col = grep('padj', colnames(res), value = T)

  res = res[order(res[,fc_col], decreasing = T),]
  genes = res[[fc_col]]
  names(genes) = res$name

  return(genes)
}



perform_GSEA = function(res, genesets){

  rg = get_ranked_genes(res)

  fgseaRes = fgsea(pathways = genesets,
                   stats    = rg,
                   minSize  = 15,
                   maxSize  = 500,
                   eps = 0)

  fgseaRes = fgseaRes[order(fgseaRes$padj),]
  return(fgseaRes)
}

## RCisTarget

download_motifRankings = function(file = c('500bp', '10kb'), destination = NULL){

  if (is.null(destination)){destination = getwd()}

  path = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/'

  if (file == '500bp'){
    fi = '/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
    path = paste0(path, fi)
  }
  else{
    fi = '/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather'
    path = paste0(path, fi)
  }

  destination = paste0(destination, fi)

  message('Starting download, might take some time... \n\n\n')
  download(path, destination, mode = 'wb')
}



load_motifData = function(rankings_path){
  motifRankings = importRankings(rankings_path)
}



perform_motif_analysis = function(res, motifRankings, motifAnnotations, name = NULL){
  #uses the RcisTarget package

  genes = res[res$significant, 'name']
  genes = list(genes)
  if (is.null(name)){names(genes) = 'genes'}
  else {names(genes) = name}

  motifs = cisTarget(genes, motifRankings, motifAnnotations)
  return(motifs)
}



show_enriched_motifs = function(motif_analysis, top_n = 10){

  resultsSubset = motif_analysis[1:top_n,]
  showLogo(resultsSubset)
}



show_TF_logo = function(motifAnnotations, TF){

  motifAnnotations = as.data.frame(motifAnnotations)
  showLogo(motifAnnotations[motifAnnotations$directAnnotation & motifAnnotations$TF == TF,])
}

