library(data.table)

setwd("~/GSE254991/database/TCGA_BRCA")

url_bulk_rna <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap/HiSeqV2.gz"
url_phenotype <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival/BRCA_survival.txt"

dest_file_bulk_rna <- "TCGA.BRCA.bulkRNA.gz"
dest_file_phenotype <- "TCGA.BRCA.survival.txt"

if (!file.exists(dest_file_bulk_rna)) {
  message("File missing. Downloading file...")
  download.file(url_bulk_rna, destfile = dest_file_bulk_rna, method = "curl")
  message("Download complete.")
} else {
  message("File '", dest_file_bulk_rna, "' already exists. Skipping download.")
}

if (!file.exists(dest_file_phenotype)) {
  message("File missing. Downloading file...")
  download.file(url_phenotype, destfile = dest_file_phenotype, method = "curl")
  message("Download complete.")
} else {
  message("File '", dest_file_phenotype, "' already exists. Skipping download.")
}

####################################################################################################
##### Load Phenotype Data #####

# Read the data
phenotype <- fread("TCGA.BRCA.survival.txt")

# Keep only 'sample' and 'OS.time' columns, and rename 'OS.time' to 'time'
phenotype <- phenotype[, .(sample, time = OS.time)]

# Create a new column 'status': assign 0 if 'sample' ends with '-11' (normal tissue), else 1 (tumor)
phenotype[, status := ifelse(grepl("-11$", sample), 0, 1)]

# Convert the data.table 'phenotype' into a standard data frame
phenotype <- as.data.frame(phenotype)

# Set the 'sample' column as row names
rownames(phenotype) <- phenotype$sample

# Remove the redundant 'sample' column from the data frame
phenotype$sample <- NULL

# head(phenotype)
# View(phenotype)
dim(phenotype)

phenotype <- phenotype[!is.na(phenotype$time), ]
sum(is.na(phenotype$time))

####################################################################################################
##### Load Bulk RNA Data #####

# Read the data
bulk_rna <- fread("TCGA.BRCA.bulkRNA.gz")

# Convert the data.table 'bulk_rna' into a standard data frame
bulk_rna <- as.data.frame(bulk_rna)

# Set the 'sample' column as row names
rownames(bulk_rna) <- bulk_rna$sample

# Remove the redundant 'sample' column from the data frame
bulk_rna$sample <- NULL

# head(bulk_rna)
# View(bulk_rna)
dim(bulk_rna)

####################################################################################################
##### Find common samples #####

common_samples <- intersect(rownames(phenotype), colnames(bulk_rna))

# Subset phenotype and bulk RNA
phenotype <- phenotype[common_samples, ]
bulk_dataset <- bulk_rna[, common_samples]

# Check
all(rownames(phenotype) == colnames(bulk_dataset))
dim(bulk_dataset)

saveRDS(
  list(
    phenotype = phenotype,
    bulk_dataset = bulk_dataset
  ),
  file = "TCGA_BRCA.rds"
)




