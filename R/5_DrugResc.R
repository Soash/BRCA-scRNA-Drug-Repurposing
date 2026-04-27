library(Seurat)
library(DrugReSC)
library(dplyr)
library(qs2)

# ==================================================================

rm(list = ls())
gc(); gc()

setwd("~/GSE254991")

result_dir <- "result/5_drug_discovery"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ------------------------------------------------------------------
### PHASE 5: DRUG DISCOVERY (DrugReSC)
### ------------------------------------------------------------------
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "STARTING PHASE 5: DrugReSC ANALYSIS\n")

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 1: Loading scPAS-scored data...\n")
combined <- qs_read("result/4_scPAS/4_annotation.qs2")


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 2: Preparing input for DrugReSC...\n")
expr_matrix <- Seurat::GetAssayData(combined, assay = "imputation", layer = "data")
PACSI_result_data <- data.frame(
  p.value = combined@meta.data$scPAS_Pvalue,
  FDR = combined@meta.data$scPAS_FDR,
  cell_phenotype_labels = ifelse(combined@meta.data$scPAS == "scPAS+", 1, 0)
)
rownames(PACSI_result_data) <- rownames(combined@meta.data)
PACSI_result <- list(result = PACSI_result_data)

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Diagnostic: Phenotype Distribution\n")
print(table(PACSI_result$result$cell_phenotype_labels))

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Memory cleanup: Removing Seurat object...\n")
rm(combined)
gc(); gc()


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 3: Loading drug signatures and running DrugReSC...\n")
drug_signature_list <- readRDS("database/DrugReSC/drug_signature_list.rds")

candidate_drugs <- DrugReSC::DrugReSC(
  expr_matrix,
  drug_signature_list,
  PACSI_result,
  D2C_transformation = "AUCell",
  ncore = 12
)

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 4: Saving drug discovery reports...\n")
write.csv(candidate_drugs, file = file.path(result_dir, "candidate_drugs.csv"), row.names = TRUE)


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "PHASE 5 COMPLETE! Results in:", result_dir, "\n")





sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()





