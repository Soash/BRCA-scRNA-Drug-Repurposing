apply_dynamic_filter <- function(seurat_obj, n_mads = 3) {
  
  # --- 1. Calculate Thresholds ---
  
  # percent.mt
  med_mt <- median(seurat_obj$percent.mt)
  mad_mt <- mad(seurat_obj$percent.mt)
  if(mad_mt == 0) mad_mt <- 1e-6 
  max_mt <- min(med_mt + (n_mads * mad_mt), 20)
  
  # nFeature_RNA (Log10 scale)
  log_feat <- log10(seurat_obj$nFeature_RNA)
  med_log_feat <- median(log_feat)
  mad_log_feat <- mad(log_feat)
  if(mad_log_feat == 0) mad_log_feat <- 1e-6
  min_feat <- max(10^(med_log_feat - (n_mads * mad_log_feat)), 200) 
  max_feat <- 10^(med_log_feat + (n_mads * mad_log_feat))
  
  # nCount_RNA (Log10 scale)
  log_count <- log10(seurat_obj$nCount_RNA)
  med_log_count <- median(log_count)
  mad_log_count <- mad(log_count)
  if(mad_log_count == 0) mad_log_count <- 1e-6
  min_count <- max(10^(med_log_count - (n_mads * mad_log_count)), 500) 
  max_count <- 10^(med_log_count + (n_mads * mad_log_count))
  
  # --- 2. Subset the object ---
  filtered_obj <- subset(seurat_obj, 
                         subset = nFeature_RNA > min_feat & 
                           nFeature_RNA < max_feat & 
                           nCount_RNA > min_count &
                           nCount_RNA < max_count &
                           percent.mt < max_mt)
  
  # --- 3. Package the Thresholds for Logging ---
  thresholds_used <- data.frame(
    mt_limit = round(max_mt, 2),
    feat_min = round(min_feat, 0),
    feat_max = round(max_feat, 0),
    count_min = round(min_count, 0),
    count_max = round(max_count, 0)
  )
  
  # Return as a list
  return(list(obj = filtered_obj, stats = thresholds_used))
}

