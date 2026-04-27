library(Seurat)
library(scPAS)
library(dplyr)
library(ggplot2)
library(tidyr)
library(qs2)

# ==================================================================

rm(list = ls())
gc(); gc()

setwd("~/GSE254991")

result_dir <- "result/4_scPAS"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

plot_dir <- file.path(result_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### PHASE 4: scPAS RISK PHENOTYPE SCORING
### ==================================================================

cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING PHASE 4: scPAS RISK SCORING\n")

combined <- qs_read("result/3_annotation/3_annotation.qs2")
TCGA_BRCA <- readRDS("database/TCGA_BRCA/TCGA_BRCA.rds")


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "scPAS start\n")
combined <- scPAS(bulk_dataset = TCGA_BRCA$bulk_dataset,
                  phenotype = TCGA_BRCA$phenotype,
                  sc_dataset = combined,
                  assay = 'SCT',
                  imputation = TRUE,
                  nfeature = 3000,
                  alpha = 0.01,
                  network_class = 'SC',
                  family = 'cox')
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "scPAS end\n")


cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Generating Visualizations...\n")

raw_pdf <- file.path(plot_dir, "Phase4_scPAS_Results.pdf")
pdf(raw_pdf, width = 12, height = 6)

# UMAP of Risk
print(DimPlot(combined, group.by = 'scPAS',
              cols = c("0" = 'grey80', "scPAS+" = 'indianred1', "scPAS-" = 'royalblue')) + 
        ggtitle("scPAS Risk Phenotype Distribution"))

# Bar Chart of Risk by Cell Type
pas_df <- as.data.frame(table(combined$SingleR_labels, combined$scPAS))
colnames(pas_df) <- c("CellType", "scPAS_Status", "Count")
pas_summary <- pas_df %>%
  group_by(CellType) %>%
  mutate(Proportion = Count / sum(Count), TotalCells = sum(Count)) %>%
  filter(TotalCells > 50) 

p_bar <- ggplot(pas_summary, aes(x = reorder(CellType, Count), y = Proportion, fill = scPAS_Status)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +
  scale_fill_manual(values = c("0" = "grey80", "scPAS+" = "indianred1", "scPAS-" = "royalblue")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "Proportion of Cells", fill = "scPAS Status",
       title = "Proportion of Risk Phenotypes by Cell Type")
print(p_bar)

dev.off()








cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 4: Compiling Detailed Risk Reports...\n")

# 1. Full Metadata Extraction
metadata_full <- combined@meta.data %>%
  select(orig.ident, TissueType, SingleR_labels, scPAS)
write.csv(metadata_full, file = file.path(result_dir, "Full_Cell_Risk_Metadata.csv"), row.names = TRUE)

# 2. Detailed Cell Type Summary
celltype_report <- metadata_full %>%
  group_by(SingleR_labels, scPAS) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(SingleR_labels) %>%
  mutate(percent = (n / sum(n)) * 100) %>%
  # Reshape to make it a one-row-per-celltype summary
  pivot_wider(names_from = scPAS, 
              values_from = c(n, percent), 
              values_fill = 0)

write.csv(celltype_report, file = file.path(result_dir, "Summary_Risk_by_CellType.csv"), row.names = FALSE)

# 3. Patient/Sample High-Risk Ranking
sample_report <- metadata_full %>%
  group_by(orig.ident, scPAS) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percent_of_sample = (n / sum(n)) * 100) %>%
  filter(scPAS == "scPAS+") %>%
  select(Sample = orig.ident, HighRisk_Count = n, HighRisk_Percent = percent_of_sample) %>%
  arrange(desc(HighRisk_Percent))

write.csv(sample_report, file = file.path(result_dir, "Summary_Risk_by_Sample.csv"), row.names = FALSE)

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Reports successfully written to:", result_dir, "\n")




qs_save(combined, file = file.path(result_dir, "4_annotation.qs2"), nthreads = 12)


sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()




