rm(list = ls())
gc(); gc()

# ==================================================================
# 1. Environment Setup
# ==================================================================
library(Seurat)
library(Asgard)
library(dplyr)
library(qs2)

setwd("~/GSE254991")

result_dir <- "result/6_Asgard"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### ASGARD PIPELINE: ORTHOGONAL DRUG VALIDATION
### ==================================================================
cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING ASGARD PIPELINE\n")

# --- Step 1: Load Data & Prepare Assay ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 1: Loading scPAS-scored Seurat object...\n")
combined <- qs_read("result/4_scPAS/4_annotation.qs2") # Loading the final object from Phase 4

# ASGARD requires standard RNA normalization for FindMarkers to work correctly
DefaultAssay(combined) <- "RNA"

cat(format(Sys.time(), "[%H:%M:%S]"), "Joining layers for Differential Expression...\n")
combined <- JoinLayers(combined)
combined <- NormalizeData(combined, verbose = FALSE)

# --- Step 2: Calculate Cell-Type Specific DEGs (Tumor vs Normal) ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 2: Calculating DEGs per Cell Type...\n")

Gene.list <- list()
C_names <- NULL
cell_types <- unique(combined$SingleR_labels)

for(i in cell_types){
  Idents(combined) <- "SingleR_labels"
  # Subset to the specific cell type
  c_cells <- subset(combined, idents = i)
  
  # Check if this cell type has enough cells in both conditions
  counts <- table(c_cells$TissueType)
  if("Tumor" %in% names(counts) && "Normal" %in% names(counts) && counts["Tumor"] > 10 && counts["Normal"] > 10){
    
    cat("Processing:", i, "... ")
    Idents(c_cells) <- "TissueType"
    
    # Calculate DEGs (Tumor vs Normal)
    C_data <- FindMarkers(c_cells, ident.1 = "Tumor", ident.2 = "Normal", verbose = FALSE)
    
    # Format for ASGARD requirements
    Gene.list[[i]] <- data.frame(
      row.names = row.names(C_data),
      score = C_data$avg_log2FC,
      adj.P.Val = C_data$p_val_adj,
      P.Value = C_data$p_val
    )
    C_names <- c(C_names, i)
    cat("Done.\n")
  } else {
    cat("Skipping", i, "- Insufficient cells in one or both groups.\n")
  }
}
names(Gene.list) <- C_names

# --- Step 3: Drug Repurposing (Mono-Drug) ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 3: Loading Breast Reference & Running Repurposing...\n")

my_gene_info <- read.table("database/L1000/DrugReference/breast_gene_info.txt", sep="\t", header = TRUE, quote = "")
my_drug_info <- read.table("database/L1000/DrugReference/breast_drug_info.txt", sep="\t", header = TRUE, quote = "")

drug.ref.profiles <- GetDrugRef(
  drug.response.path = 'database/L1000/DrugReference/breast_rankMatrix.txt',
  probe.to.genes = my_gene_info, 
  drug.info = my_drug_info
)

Drug.ident.res <- GetDrug(
  gene.data = Gene.list, 
  drug.ref.profiles = drug.ref.profiles, 
  repurposing.unit = "drug", 
  connectivity = "negative", 
  drug.type = "FDA"
)

# 35 min
qs_save(Drug.ident.res, file = file.path(result_dir, "Drug.ident.res.qs2"), nthreads = 12)
# Drug.ident.res <- qs_read(file.path(result_dir, "Drug.ident.res.qs2"))

# --- Step 4: Final Scoring & Fixed Export ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 4: Calculating Therapeutic Scores & Saving CSV...\n")

cell_metadata <- combined@meta.data
cell_metadata$cluster <- combined$SingleR_labels
cell_metadata$sample <- combined$orig.ident 
Case_samples <- unique(cell_metadata$sample[cell_metadata$TissueType == "Tumor"])

Drug.score <- DrugScore(
  cell_metadata, 
  cluster_degs = Gene.list, 
  cluster_drugs = Drug.ident.res, 
  tissue = "breast", 
  case = Case_samples, 
  gse92742_gctx_path = "database/L1000/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", 
  gse70138_gctx_path = "database/L1000/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
)

# Convert row names to a permanent column and filter
Final.drugs <- Drug.score %>%
  mutate(Drug_Name = rownames(.)) %>%
  filter(FDR < 0.05) %>%
  arrange(desc(Drug.therapeutic.score)) %>%
  dplyr::select(Drug_Name, Drug.therapeutic.score, P.value, FDR) # Added dplyr:: here

# Save the clean CSV
write.csv(Final.drugs, file.path(result_dir, "ASGARD_Final_Mono_Drugs.csv"), row.names = FALSE)

cat(format(Sys.time(), "[%H:%M:%S]"), "ASGARD PIPELINE COMPLETE!\n")

# --- Teardown ---
sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()







