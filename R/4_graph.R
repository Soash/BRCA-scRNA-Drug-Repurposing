library(Seurat)
library(dplyr)
library(ggplot2)
library(qs2)

# 1. Load the already-scored Phase 4 object
combined <- qs_read("result/4_scPAS/4_annotation.qs2")
plot_dir <- "result/4_scPAS/plots"

# 2. Generate and save the UMAP as a high-res PNG (Figure 2A)
p_umap <- DimPlot(combined, group.by = 'scPAS',
                  cols = c("0" = 'grey80', "scPAS+" = 'indianred1', "scPAS-" = 'royalblue')) + 
  ggtitle("scPAS Risk Phenotype Distribution")

ggsave(file.path(plot_dir, "Figure2A_Risk_UMAP.png"), plot = p_umap, width = 8, height = 6, dpi = 300)

# 3. Generate and save the Bar Chart as a high-res PNG (Figure 2B)
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

ggsave(file.path(plot_dir, "Figure2B_Risk_BarChart.png"), plot = p_bar, width = 8, height = 8, dpi = 300)

print("High-resolution PNGs saved to result/4_scPAS/plots/")