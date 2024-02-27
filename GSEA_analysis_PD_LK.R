### GSEA analysis with fgsea package from Bioconductor
# Requires the experimental design in csv file
# Requires the pg.txt file exported from pd.result file
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(enrichR)

# Load data -> Import the data using the "Import Dataset" from the environment
# Rename the raw data analysis PD output txt file to "data"
# Rename the  ExperimentalDesign.csv file to "ExperimentalDesign"


# Read all data 

ExperimentalDesign$replicate <- as.numeric(ExperimentalDesign$replicate )

# Convert columns that match the regex (basically the columns that contain maxLFQ intensities or equivalent) to numeric
columns_to_convert <- grep("^Abundances..Grouped...F\\d", names(data), value = TRUE)
data[columns_to_convert] <- lapply(data[columns_to_convert], as.numeric)

# Function to extract Gene Name from Description
extract_gene_name <- function(description) {
  # Find the position of "GN=" in the Description
  gn_position <- regexpr("GN=", description)
  
  if (gn_position != -1) {
    # Extract the substring after "GN=" until the next space
    gene_name <- sub(".*GN=([^\\s]+).*", "\\1", description)
    gene_name <- sub(" PE=\\d+ SV=\\d+", "", gene_name)
    return(gene_name)
  } else {
    # If "GN=" is not found, return NA or handle accordingly
    return(NA)
  }
}

# Apply the function to the Description column and create a new column 'Gene.names'
data$Gene.names <- sapply(data$Description, extract_gene_name)

# Rename the column "Accession" to "Protein.IDs"
colnames(data)[colnames(data) == "Accession"] <- "Protein.IDs"


# filter the data if necessary (for MQ)
# data <- filter(data, Reverse != "+", Potential.contaminant != "+")


#Verify any duplicate Gene.names
data$Genes.names %>% duplicated() %>% any()

# Make a table of duplicated gene names
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
# Verification to check duplicity. Aim: FALSE
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("^Abundances..Grouped...F\\d", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se(data_unique, LFQ_columns, ExperimentalDesign)
# Generate a SummarizedExperiment object by parsing condition information from the column names
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)

# Filter, normalize and impute missing values
data_filt <- filter_missval(data_se, thr = 0)
data_norm <- normalize_vsn(data_filt)
data_imputed <- impute(data_norm, fun = "MinProb", q = 0.01)

# Test for differentially expressed proteins
diff <- test_diff(data_imputed, type='all') 
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)


# DO the GSEA analysis 
# you can change the database
gsea_kegg <-test_gsea(dep, databases = "GO_Cellular_Component_2017b", contrasts = TRUE) # Logical Whether or not to perform the gene set enrichment analysis independently for the different contrasts.

### Visualization

ggplot_object <- plot_gsea(gsea_kegg, number = 10, alpha = 0.05, contrasts = NULL,
                           databases = NULL, nrow = 1, term_size = 8)  
print(ggplot_object)
