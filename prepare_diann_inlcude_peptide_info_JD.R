library(reshape2)
library(diann)
library(dplyr)

pg = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/V63_Paulien/report.pg_matrix_T_cells.tsv', stringsAsFactors = F, sep = '\t') # = pg.matrix
rep = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/V63_Paulien/report_T_cells.tsv', stringsAsFactors = F, sep = '\t') # = report
contaminants = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/contaminants.txt', sep = '\t')$x #common contaminants file always in the same location

#pg = data2 #pg matrix
#rep = data1 #report


getPeptideNo = function(rep){

  df = dcast(rep, Precursor.Id + Protein.Group + Precursor.Id + Proteotypic ~ File.Name, value.var = 'Precursor.Normalised')
  summed_int = aggregate(as.matrix(df[,5:10]) ~ Precursor.Id, data = df, FUN = sum, na.rm=TRUE, na.action=NULL)
  summed_int$Protein.Group = df[match(summed_int$Precursor.Id, df$Precursor.Id), 'Protein.Group']
  summed_int$Proteotypic = df[match(summed_int$Precursor.Id, df$Precursor.Id), 'Proteotypic']

  l_total = split(summed_int, summed_int$Protein.Group)
  v_total = sapply(l_total, function(x){
    ints = x[,2:7]
    pep = colSums(!(ints == 0))
    pep = c(pep, total = nrow(ints))
  })

  return(t(v_total))
}

prepare_pg = function(pg, report, contaminants){

  df = report
  protein.groups = diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")

  protein.groups = (function(protein.groups, report){

    id_map = report[!duplicated(report$Protein.Group),]
    id_map = id_map[,3:6]

    protein.groups = merge(id_map, protein.groups, by.x = 'Protein.Group', by.y = 0)
    return(protein.groups)
  })(protein.groups, report)

  protein.groups$Potential.contaminant = ifelse(protein.groups$Protein.Group %in% contaminants, '+', '')

  pep = getPeptideNo(report)
  protein.groups = merge(protein.groups, pep, by.x = 'Protein.Group',  by.y = 0, all.x = T, all.y = F)

  return(protein.groups)
}


#the final column of this df contains information about the number of unique peptides found across all samples for a particular protein group
pg = prepare_pg(pg, rep, contaminants)

#filter the dataframe to keep proteins that have at least 2 unique peptides 
# adjust if you want to be more stringent
data <- pg %>%
  filter(total >= 2)

#drop the columns that you don't need anymore (so you can feed it to dep package)

data <- data %>%
  select(1:16) #change the number each time by calculating: 4 standards columns + number of samples you have 

# Save the combined dataframe as a tsv file in the specified directory ->always change the file name 
data_folder <- "R:/Group Vermeulen/Lila/Mass_spec_results/V63_Paulien"

write.table(data, file = paste0(data_folder, "/report.pg_peptides_filtered.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
