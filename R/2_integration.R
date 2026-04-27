library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(qs2)

# ==================================================================

rm(list = ls())
gc(); gc()

setwd("~/GSE254991")

result_dir <- "result/2_integration"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

plot_dir <- file.path(result_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### PHASE 2: MERGING & HARMONY INTEGRATION
### ==================================================================

cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING PHASE 2: INTEGRATION\n")


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 1: Load & Merge\n")
seurat_list <- qs_read("result/1_qc/1_QC.qs2")
combined <- merge(seurat_list[[1]], y = seurat_list[-1])


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 2: Basic Stats\n")
print(combined)
count_genes <- sum(rowSums(combined[["RNA"]]$counts) > 0)
count_cells <- ncol(combined)
cat("Genes:", count_genes, "| Cells:", count_cells, "\n")
cat("Total Samples Integrated: ", length(unique(combined$orig.ident)), "\n")
rm(seurat_list); gc(); gc()

options(future.globals.maxSize = 20 * 1024^3)
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 3: Running SCTransform...\n")
combined <- SCTransform(combined)


cat("\n\n", format(Sys.time(), "[%H:%M:%S]"), "Step 4: Running PCA...\n")
combined <- RunPCA(combined)


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 5: Running Harmony Integration...\n")
combined <- IntegrateLayers(
  object = combined,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca", 
  new.reduction = "integrated.dr",
  verbose = TRUE
)


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 6: Running UMAP and Clustering...\n")
combined <- FindNeighbors(combined, reduction = "integrated.dr", dims = 1:30)
combined <- FindClusters(combined)
combined <- RunUMAP(combined, reduction = "integrated.dr", dims = 1:30)


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 7: Saving Phase 2 Plots...\n")
raw_pdf <- file.path(plot_dir, "Phase2_Integration_UMAPs.pdf")
pdf(raw_pdf, width = 12, height = 6)
print(DimPlot(combined, group.by = "orig.ident") + ggtitle("Integrated UMAP by Patient"))
print(DimPlot(combined, group.by = "TissueType") + ggtitle("Integrated UMAP by Tissue Type"))
print(DimPlot(combined, label = TRUE) + ggtitle("Integrated UMAP by Cluster"))
dev.off()


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 8: Saving integrated qs2...\n")
qs_save(combined, file = file.path(result_dir, "2_integration.qs2"), nthreads = 12)
cat(format(Sys.time(), "[%H:%M:%S]"), "PHASE 2 COMPLETE!\n")


sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()

