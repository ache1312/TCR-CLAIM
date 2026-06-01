read_tcrbio_table <- function(path) {
  if (!file.exists(path)) {
    stop("Table does not exist: ", path, call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tcrbio_table <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE, na = "")
  invisible(path)
}

write_tcrbio_tables <- function(cell_tcr, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_tcrbio_table(cell_tcr, file.path(out_dir, "cell_tcr_table.csv"))
}

as_tcrbio_cell_table <- function(cell_metadata_with_tcr) {
  df <- as.data.frame(cell_metadata_with_tcr, stringsAsFactors = FALSE)
  rename_if_present <- function(old, new) {
    if (old %in% names(df) && !new %in% names(df)) {
      names(df)[names(df) == old] <<- new
    }
  }
  rename_if_present("tcr_cell_barcode", "cell_barcode")
  required <- c(
    "cell_barcode", "dataset_id", "donor_id", "sample_id", "tissue_type",
    "cell_class", "ct_strict", "ct_vgene"
  )
  for (column in required) {
    if (!column %in% names(df)) {
      df[[column]] <- NA
    }
  }
  df
}
