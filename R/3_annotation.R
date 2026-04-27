library(Seurat)
library(SingleR)
library(BiocParallel)
library(ggplot2)
library(dplyr)
library(qs2)

# ==================================================================

rm(list = ls())
gc(); gc()

setwd("~/GSE254991")

result_dir <- "result/3_annotation"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

plot_dir <- file.path(result_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### PHASE 3: AUTOMATED CELL TYPE ANNOTATION (SingleR)
### ==================================================================


cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING PHASE 3: ANNOTATION\n")


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 1: Loading integrated qs2 file...\n")
combined <- qs_read("result/2_integration/2_integration.qs2")


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 2: Loading THBCA reference...\n")
THBCA <- readRDS("database/THBCA/THBCA_max500_per_celltype.rds")


overlap_genes <- intersect(rownames(combined), rownames(THBCA$ref))
cat(format(Sys.time(), "[%H:%M:%S]"), "Gene overlap found:", length(overlap_genes), "genes.\n")


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 3: Running SingleR (Parallel on 12 cores)...\n")
test_data <- GetAssayData(combined, assay = "SCT", layer = "data")
pred_cell <- SingleR(
  test = test_data,
  ref = THBCA$ref,
  labels = THBCA$cell_type,
  de.method = "wilcox",
  BPPARAM = MulticoreParam(workers = 12, progressbar = TRUE)
)


combined$SingleR_labels <- pred_cell$labels # Assign labels back to Seurat metadata
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "SingleR finished in")


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 8: Saving annotation qs2...\n")
qs_save(combined, file = file.path(result_dir, "3_annotation.qs2"), nthreads = 12)
cat(format(Sys.time(), "[%H:%M:%S]"), "PHASE 3 COMPLETE!\n")


sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()





