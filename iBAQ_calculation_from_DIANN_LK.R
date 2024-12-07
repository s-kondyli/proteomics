# Load required libraries
library(Peptides)
library(stringr)
library(dplyr)
library(data.table)
# requires as input pg & pr matrices

# Extract and save as tsv the protein group ids from the pg matrix so you can run the "uniprot batch retrieval" online
pg = read.delim('R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all/all_filtered_2_unique_peptides/report.pg_peptides_filtered_all_PBMCs_DIANN_1.9.2 - Copy.tsv', stringsAsFactors = F, sep = '\t') # = pg.matrix
# Clean up the 'Protein.Group' column to keep only the part before the first ';'
pg <- pg %>%
  mutate(Protein.Group = sub(";.*", "", Protein.Group))

# Define the output path (replace this with your desired file path)
data_folder <- "R:/Group Vermeulen/Lila/Mass_spec_results/PBMCs_all/all_filtered_2_unique_peptides"

# Save only the 'Protein.Group' column as a TSV file -> you will use it for "uniprot batch retrieval" online
write.table(pg$Protein.Group, file = paste0(data_folder, "/Protein_IDs.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Perform uniprot batch retrieval online and download the file in tsv format -> make sure to check the "Sequence" column under tab "Sequences"
# Read the tsv file you downloaded from Uniprot as idmapping and make sure to check the "Sequence" column is present
# With this process you get the protein sequences for all your proteins identified in your dataset
# Necessary for in silico digestion tailored to your dataset -> instead of using a FASTA file which would also 
#contain proteins that are not present in your dataset and you would need to filter for only the proteins present in your dataset 

# Function to perform in silico digestion with trypsin and apply filtering criteria
digest_protein <- function(protein_sequence) {
  # Define trypsin cleavage rule (cleaves after K and R, but not if followed by P)
  cleavage_sites <- gregexpr("(?<=[KR])(?!P)", protein_sequence, perl = TRUE)[[1]] + 1
  
  # If no cleavage sites, return 0 peptides
  if (length(cleavage_sites) == 0) {
    return(0)
  }
  
  # Generate peptides
  peptides <- c()
  start_pos <- 1
  for (i in cleavage_sites) {
    peptides <- c(peptides, substr(protein_sequence, start_pos, i - 1))
    start_pos <- i
  }
  peptides <- c(peptides, substr(protein_sequence, start_pos, nchar(protein_sequence)))  # Add final peptide
  
  # Filter peptides by length (6-30) 
  valid_peptides <- peptides[nchar(peptides) >= 6 & nchar(peptides) <= 30]
  
  # Use the mw() function from Peptides package to calculate peptide masses
  peptide_masses <- sapply(valid_peptides, function(pep) {
    # Calculate the mass using mw() from the Peptides package
    mw(pep, monoisotopic = FALSE, avgScale = "expasy", label = "none", aaShift = NULL)
  })
  
  # Filter peptides by mass (500-5000 Da)
  valid_peptides <- valid_peptides[peptide_masses >= 500 & peptide_masses <= 5000]
  
  # Return the number of valid peptides
  return(length(valid_peptides))
}

# Applying this function to a dataframe with a column of protein sequences
apply_peptide_count <- function(df, sequence_column) {
  df$observable_peptide_count <- sapply(df[[sequence_column]], digest_protein)
  return(df)
}

#Apply the function which creates a new column with the theoretically observable peptides 
# first argument of the function expects the dataframe to be processed
# second argument of the function expects the column name of the dataframe that contains the protein sequences

idmapping <- apply_peptide_count(idmapping, "Sequence")

# Calculate the mean of precursor normalized intensity for each peptide accross all samples
# Input: report.pr_matrix
# Convert your dataframe to a data.table for faster processing because of a lot of data
setDT(report.pr_matrix)

# Perform the grouping and averaging using data.table
df_avg_intensity <- report.pr_matrix[, 
                                     lapply(.SD, mean, na.rm = TRUE),  # Compute mean for each column in .SD (Subset of Data)
                                     by = .(Protein.Group, Stripped.Sequence, Genes),  # Group by these columns
                                     .SDcols = 11:ncol(report.pr_matrix)  # Select only the sample columns starting from the 11th column
                                                                          # because always the pr_matrix samples start from the 11th column
]

# Summarize peptide intensity for each Protein Group in each sample by summing the mean of the precursor normalized intensities
df_avg_protein_intensity <- df_avg_intensity[, 
                                             lapply(.SD, sum, na.rm = TRUE),  # Sum the precursor normalized intensities for each sample
                                             by = .(Protein.Group, Genes),  # Group by Protein.Group (you can add other grouping variables like Genes if needed)
                                             .SDcols = 4:ncol(df_avg_intensity)  # Only sum over the sample columns 
                                                                                 # starting with the 4rth column because the 
                                                                                 # input is the df_avg_intensity data table
]

#Merge the information between idmapping & df_avg_protein_intensity

merged_data <- df_avg_protein_intensity %>%
  left_join(idmapping, by = c("Protein.Group" = "Entry"))

# Calculate iBAQ for each protein in each sample
# Select sample columns -> change each time
sample_columns <- grep("^F..LK", colnames(merged_data), value = TRUE)

# Apply the iBAQ formula to each sample
merged_data <- merged_data %>%
  mutate(across(all_of(sample_columns), 
                ~ . / observable_peptide_count, 
                .names = "{.col}_iBAQ"))

# Save your resultrs
write.table(merged_data, file = paste0(data_folder, "/iBAQ_info_PBMCs_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
