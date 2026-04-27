library(DrugReSC)

setwd("~/GSE254991/database/DrugReSC")

source("get_drug_iLINCS_sig_safe.R")

drug_iLINCS_sig_id <- c(
  "LINCSCP_178870", "LINCSCP_179309", "LINCSCP_2385", "LINCSCP_2554", 
  "LINCSCP_26", "LINCSCP_34", "LINCSCP_53", "LINCSCP_53344", 
  "LINCSCP_53471", "LINCSCP_69009"
)

iLINCS_signature_list <- get_drug_iLINCS_sig_safe(drug_iLINCS_sig_id)

successful_downloads <- sum(!sapply(iLINCS_signature_list, is.null))
cat("Successfully downloaded", successful_downloads, "out of", length(drug_iLINCS_sig_id), "signatures.\n")

drug_signature_list <- DrugReSC::get_drug_sig(iLINCS_signature_list)

saveRDS(drug_signature_list, file = "drug_signature_list.rds")








