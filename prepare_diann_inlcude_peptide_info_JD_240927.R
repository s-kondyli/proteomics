
### WARNING WARNING WARNING###
### HAS NOT BEEN TESTED EXTENSIVELY


library(reshape2)
library(diann)
library(dplyr)
library(diann)

pg = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all/report.pg_matrix.tsv', stringsAsFactors = F, sep = '\t') # = pg.matrix
rep = read.delim("R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all/report.tsv", stringsAsFactors = F, sep = '\t') # = report
contaminants = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/contaminants.txt', sep = '\t')$x #common contaminants file always in the same location
#alternative for opening the report file fast
#rep =  diann_load("R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all/report.tsv")

getPeptideNo = function(rep, id_column = 'Protein.Group', peptide_column = 'Precursor.Id', quantity_column = 'Precursor.Normalised'){
  
  
  form = paste(peptide_column, id_column, 'Proteotypic', sep = '+')
  form = as.formula(paste(form, '~ File.Name'))
  
  df = dcast(rep, form, value.var = quantity_column, fun.aggregate = sum)
  summed_int = df
  summed_int$Protein.Group = df[match(summed_int[[peptide_column]], df[[peptide_column]]), id_column]
  summed_int$Proteotypic = df[match(summed_int[[peptide_column]], df[[peptide_column]]), 'Proteotypic']
  
  l_total = split(summed_int, summed_int[[id_column]])
  v_total = sapply(l_total, function(x){
    ints = x[,4:(ncol(x)-1)]
    pep = colSums(!(ints == 0))
    pep = c(pep, total = nrow(ints))
  })
  
  return(t(v_total))
}


prepare_pg = function(pg, report, contaminants, id_column = 'Protein.Group', quantity_column = 'Precursor.Normalised', id_header = 'Precursor.Id',
                      npep_level = 'Precursor.Id', npep_quantity_col = 'Precursor.Normalised' ){
  
  # id_column specifies which column will be used as rownames in the final data
  # quantity_column specifies which column is used to determine the LFQ intensities
  # id_header specifies the column with the peptide id information.
  # npep_level specifies the peptide id column to determine the number or peptides detected. In case of duplicates, intensities are summed
  # npep_quantity_col specifies which quantity column is used to determine peptide number. If a value is NA, the peptide is not counted.

  df = report
  protein.groups = diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header=id_column, id.header = id_header, quantity.header = quantity_column)

  protein.groups = (function(protein.groups, report){

    id_map = report[!duplicated(report$Protein.Group),]
    id_map = id_map[,3:6]

    protein.groups = merge(id_map, protein.groups, by.x = id_column, by.y = 0)
    return(protein.groups)
  })(protein.groups, report)

  protein.groups$Potential.contaminant = ifelse(protein.groups$Protein.Group %in% contaminants, '+', '')

  pep = getPeptideNo(report, id_column, npep_level, npep_quantity_col)
  protein.groups = merge(protein.groups, pep, by.x = id_column,  by.y = 0, all.x = T, all.y = F, suffixes = c('_LFQ_intensity', '_n_peptides'))

  return(protein.groups)
}


#the final column of this df contains information about the number of unique peptides found across all samples for a particular protein group
pg = prepare_pg(pg, rep, contaminants, 'Genes', id_header = 'Stripped.Sequence', quantity_column = 'Precursor.Normalised', npep_level = 'Stripped.Sequence', npep_quantity_col = 'Precursor.Normalised')
#pg2 = prepare_pg(pg, rep, contaminants, 'Genes', id_header = 'Precursor.Id', quantity_column = 'Precursor.Normalised', npep_level = 'Precursor.Id', npep_quantity_col = 'Precursor.Normalised')

#option to save pg2 matrix as well
# Save the combined dataframe as a tsv file in the specified directory ->always change the file name 
data_folder <- "R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all"

write.table(pg, file = paste0(data_folder, "/report.pg_peptide_info_all_PBMCs_DIANN_1.9.2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

#filter the dataframe to keep proteins that have at least 2 unique + razor peptides 
# adjust if you want to be more stringent
data <- pg %>%
  filter(total >= 2)

# Filter out the contaminants
data <- data %>%
      filter(Potential.contaminant != "+")

#drop the columns that you don't need anymore (so you can feed it to dep package)

data <- data %>%
  select(1:28) #change the number each time by calculating: 4 standards columns + number of samples you have 

# Save the combined dataframe as a tsv file in the specified directory ->always change the file name 
data_folder <- "R:/Group Vermeulen/Lila/Mass_spec_results/V_LK_One_tip"

write.table(data, file = paste0(data_folder, "/report.pg_peptides_filtered_all_PBMCs_DIANN_1.9.2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
