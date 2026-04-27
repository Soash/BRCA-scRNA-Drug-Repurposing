rm(list = ls())
gc(); gc()

# ==================================================================
# 1. Environment Setup
# ==================================================================
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(qs2)

setwd("~/GSE254991")

result_dir <- "result/6_pathway"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ------------------------------------------------------------------
### PHASE 6: PATHWAY ENRICHMENT ANALYSIS (Functional Genomics)
### ------------------------------------------------------------------
cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING PHASE 6: PATHWAY ENRICHMENT\n")

# --- Step 1: Load and Prepare Data ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 1: Loading scPAS-scored Seurat object...\n")
combined <- qs_read("result/4_scPAS/4_annotation.qs2")

# PrepSCTFindMarkers ensures the different sample models are reconciled
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 2: Reconciling SCT models for FindMarkers...\n")
combined <- PrepSCTFindMarkers(combined)

# --- Step 2: Find Markers (Risk Phenotype) ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 3: Calculating DEGs between scPAS+ and scPAS- cells...\n")
Idents(combined) <- "scPAS"

# We calculate markers for both up and down to get the full biological landscape
scPAS_markers <- FindMarkers(
  combined, 
  ident.1 = "scPAS+", 
  ident.2 = "scPAS-", 
  assay = "SCT", 
  recode.combined.layers = TRUE, # Seurat v5 specific
  logfc.threshold = 0.25, 
  min.pct = 0.1
)

write.csv(scPAS_markers, file.path(result_dir, "Full_DEGs_scPAS_plus_vs_minus.csv"))

# --- Step 3: Gene ID Conversion (SYMBOL to ENTREZ) ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 4: Filtering genes and converting IDs...\n")

# Focus on significant upregulated genes (High Risk drivers)
sig_up <- subset(scPAS_markers, p_val_adj < 0.05 & avg_log2FC > 0.5)
sig_gene_symbols <- rownames(sig_up)

cat("Found", length(sig_gene_symbols), "significant upregulated genes in scPAS+ cells.\n")

if(length(sig_gene_symbols) > 5) {
  gene_ids <- bitr(sig_gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
  
  # --- Step 4: Run Enrichment (KEGG & GO) ---
  cat(format(Sys.time(), "[%H:%M:%S]"), "Step 5: Running Enrichment Analyses...\n")
  
  kegg_results <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
  go_results <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
  
  # --- Step 5: Save & Plot Results ---
  cat(format(Sys.time(), "[%H:%M:%S]"), "Step 6: Saving results and generating plots...\n")
  
  write.csv(as.data.frame(kegg_results), file.path(result_dir, "KEGG_Enrichment_scPAS_plus.csv"), row.names = FALSE)
  write.csv(as.data.frame(go_results), file.path(result_dir, "GO_Enrichment_scPAS_plus.csv"), row.names = FALSE)
  
  pdf(file.path(result_dir, "Phase6_scPAS_Pathway_Plots.pdf"), width = 10, height = 8)
  
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    print(dotplot(kegg_results, showCategory = 15) + ggtitle("KEGG Pathways: Activated in scPAS+"))
  }
  
  if (!is.null(go_results) && nrow(go_results) > 0) {
    print(dotplot(go_results, showCategory = 15) + ggtitle("GO Biological Processes: Activated in scPAS+"))
  }
  
  dev.off()
  
} else {
  cat("WARNING: Too few significant genes found for enrichment analysis.\n")
}

cat(format(Sys.time(), "[%H:%M:%S]"), "PHASE 6 COMPLETE!\n")

# --- Teardown ---
sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()





