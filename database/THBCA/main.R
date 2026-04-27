### ==================================================================
### DATABASE PREP: THBCA Reference for SingleR
### ==================================================================

library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(Matrix)

setwd("/home/soash_sadat/GSE254991/database/THBCA")

dest_file <- "THBCA.h5ad"
n_Cell <- 500


cat(format(Sys.time(), "[%H:%M:%S]"), "Reading H5AD file...\n")
sce <- readH5AD(dest_file, use_hdf5 = TRUE, reader = "R")

cat(format(Sys.time(), "[%H:%M:%S]"), "Mapping 'feature_name' to Rownames...\n")
new_rownames <- make.unique(as.character(rowData(sce)$feature_name))
rownames(sce) <- new_rownames


cat(format(Sys.time(), "[%H:%M:%S]"), "Sampling", n_Cell, "cells per cell type...\n")
set.seed(1)
sampled_cells <- colData(sce) |>
  as.data.frame() |>
  tibble::rownames_to_column("cell_id") |>
  group_by(cell_type) |>
  slice_sample(n = n_Cell) |>  
  pull(cell_id)

sce_subset <- sce[, sampled_cells]


cat(format(Sys.time(), "[%H:%M:%S]"), "Realizing into RAM and naming assay 'logcounts'...\n")
assay(sce_subset, "logcounts") <- as(assay(sce_subset, "X"), "dgCMatrix") 
assay(sce_subset, "X") <- NULL

thbca_final <- list(
  ref = sce_subset,
  cell_type = as.character(sce_subset$cell_type)
)

rds_filename <- paste0("THBCA_max", n_Cell, "_per_celltype.rds")
saveRDS(thbca_final, file = rds_filename)
cat(format(Sys.time(), "[%H:%M:%S]"), "SUCCESS: Reference saved as", rds_filename, "\n")



