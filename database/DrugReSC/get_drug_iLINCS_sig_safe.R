get_drug_iLINCS_sig_safe <- function(drug_sig_ids) {
  iLINCS_signature_list <- list()
  for (i in seq_along(drug_sig_ids)) {
    ilincs_signatureId <- drug_sig_ids[i]
    message("Fetching ", ilincs_signatureId, " ...")
    req <- httr::POST("https://www.ilincs.org/api/ilincsR/downloadSignature",
                      body = list(sigID = ilincs_signatureId, display = FALSE),
                      encode = "json")
    content <- httr::content(req)
    if (length(content) == 0 || is.null(content[[1]])) {
      warning("Failed to fetch signature: ", ilincs_signatureId)
      iLINCS_signature_list[[i]] <- NULL
      next
    }
    ilincs_sessionId <- unlist(content)
    fileUrl <- paste0("https://www.ilincs.org/tmp/", ilincs_sessionId, ".xls")
    message("Downloading from ", fileUrl)
    tryCatch({
      signatureData <- read.table(fileUrl, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      iLINCS_signature_list[[i]] <- signatureData
    }, error = function(e) {
      warning("Could not read file for ", ilincs_signatureId, ": ", e$message)
      iLINCS_signature_list[[i]] <- NULL
    })
  }
  names(iLINCS_signature_list) <- drug_sig_ids
  return(iLINCS_signature_list)
}
