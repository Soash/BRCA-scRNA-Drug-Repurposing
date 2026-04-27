library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)
library(ggplot2)
library(future)
library(qs2)

# ==================================================================

rm(list = ls())
gc(); gc()

setwd("~/GSE254991")

result_dir <- "result/1_qc"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

plot_dir <- file.path(result_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### PHASE 1: SEQUENTIAL QC & DOUBLET REMOVAL
### ==================================================================

plan("sequential")
source("R/dynamic_filter.R")

# Get Sample Lists
cancer_samples <- list.dirs("data/cancer", full.names = TRUE, recursive = FALSE)
normal_samples <- list.dirs("data/normal", full.names = TRUE, recursive = FALSE)
all_samples <- c(cancer_samples, normal_samples)

seurat_list <- list()

# ==================================================================
# BEGIN SEQUENTIAL LOOP
# ==================================================================
for (i in seq_along(all_samples)) {
  
  x <- all_samples[i]
  sample_name <- basename(x)
  
  cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING SAMPLE", i, "of", length(all_samples), ":", sample_name, "\n")
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step A: Loading 10X data...\n")
  obj <- CreateSeuratObject(counts = Read10X(data.dir = x))
  obj$orig.ident <- sample_name
  obj$TissueType <- ifelse(grepl("cancer", x), "Tumor", "Normal")
  
  print(obj)
  before_genes <- sum(rowSums(obj[["RNA"]]$counts) > 0)
  before_cells <- ncol(obj)
  cat("Genes (Raw):", before_genes, "| Cells (Raw):", before_cells, "\n")
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step B: Calculating QC Metrics...\n")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- NormalizeData(obj, verbose = FALSE)
  print(obj)
  
  raw_pdf <- file.path(plot_dir, paste0("Raw_QC_", sample_name, ".pdf"))
  pdf(raw_pdf, width = 12, height = 6)
  print(
    VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
      ggtitle(paste(sample_name, "- Raw"))
  )
  dev.off()
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step C: Applying Dynamic 3-MAD Filter...\n")
  filter_results <- apply_dynamic_filter(obj, n_mads = 3)
  obj <- filter_results$obj
  params <- filter_results$stats
  
  print(obj)
  after_filter_genes <- sum(rowSums(obj[["RNA"]]$counts) > 0)
  after_filter_cells <- ncol(obj)
  cat("Genes (filter):", after_filter_genes, "| Cells (filter):", after_filter_cells, "\n")
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step D: Running scDblFinder...\n")
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce)
  obj$DoubletClass <- sce$scDblFinder.class
  obj <- subset(obj, subset = DoubletClass == "singlet")
  
  cat("\n\n")
  print(obj)
  final_genes <- sum(rowSums(obj[["RNA"]]$counts) > 0)
  final_cells <- ncol(obj)
  cat("Genes (Raw):", final_genes, "| Cells (Raw):", final_cells, "\n")
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step E: Finalizing logs...\n")
  cat(
    "Sample:", sample_name, "| Tissue:", obj$TissueType[1], "\n",
    "[Thresholds Applied]\n",
    "- Max MT%:", params$mt_limit, "\n",
    "- nFeature Range:", params$feat_min, "-", params$feat_max, "\n",
    "- nCount Range:", params$count_min, "-", params$count_max, "\n",
    "[Counts Tracking]\n",
    "- Cells:", before_cells, "(Raw) ->", final_cells, "(Final Singlets)\n",
    "- Genes:", before_genes, "(Raw) ->", final_genes, "(Active)\n",
    "- Doublets removed:", (after_filter_cells - final_cells), "\n",
    "--------------------------------------------------\n"
  )
  
  clean_pdf <- file.path(plot_dir, paste0("QC_", sample_name, ".pdf"))
  pdf(clean_pdf, width = 12, height = 6)
  print(
    VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
      ggtitle(paste(sample_name, "- Cleaned"))
  )
  dev.off()
  
  seurat_list[[sample_name]] <- obj
}

qs_save(seurat_list, file = file.path(result_dir, "1_QC.qs2"), nthreads = 12)

sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()


