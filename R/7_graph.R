library(ggplot2)
library(dplyr)
library(stringr)

# ==================================================================
# Standalone Script: Plot Enrichment from CSV
# ==================================================================

# Define your working directories (Adjust if necessary)
setwd("~/GSE254991")
result_dir <- "result/6_pathway"
plot_dir <- file.path(result_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- Plotting Function ---
# This function recreates the clusterProfiler dotplot perfectly from a CSV
plot_from_csv <- function(csv_path, plot_title) {
  
  # Read the CSV and take the top 15 most significant pathways
  df <- read.csv(csv_path) %>%
    head(15) %>%
    mutate(
      # Convert the "Count/Total" GeneRatio string into a numeric decimal for plotting
      GeneRatio_num = sapply(strsplit(as.character(GeneRatio), "/"), 
                             function(x) as.numeric(x[1]) / as.numeric(x[2])),
      # Wrap long pathway names so they fit nicely on the Y-axis
      Description = stringr::str_wrap(Description, width = 45),
      # Lock the order so the most significant pathways stay at the top
      Description = factor(Description, levels = rev(unique(Description)))
    )
  
  # Build the ggplot
  p <- ggplot(df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    # clusterProfiler default uses red for high significance (low p-value) and blue for lower
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    labs(title = plot_title, 
         x = "GeneRatio", 
         y = NULL, 
         size = "Count", 
         color = "p.adjust") +
    theme(
      axis.text.y = element_text(size = 11, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  return(p)
}

# --- Generate the Plots ---
cat("Reading CSVs and generating plots...\n")

p_kegg <- plot_from_csv(file.path(result_dir, "KEGG_Enrichment_scPAS_plus.csv"), 
                        "KEGG Pathways: Activated in scPAS+ cells")

p_go <- plot_from_csv(file.path(result_dir, "GO_Enrichment_scPAS_plus.csv"), 
                      "GO Biological Processes: Activated in scPAS+ cells")

# --- Save Outputs (Height = 10) ---
cat("Saving high-res outputs (Height = 10)...\n")

# 1. Save as a combined PDF
pdf(file.path(plot_dir, "Figure3_scPAS_Pathway_Plots.pdf"), width = 10, height = 10)
print(p_kegg)
print(p_go)
dev.off()

# 2. Save as individual high-res PNGs for your Thesis/Google Doc
ggsave(file.path(plot_dir, "Figure3A_KEGG_Dotplot.png"), plot = p_kegg, width = 10, height = 10, dpi = 300)
ggsave(file.path(plot_dir, "Figure3B_GO_Dotplot.png"), plot = p_go, width = 10, height = 10, dpi = 300)

cat("Done! Check the 'plots' folder for your Figure 3 images.\n")