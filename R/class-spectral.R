#' Create a RUFFS spectra object
#' 
#' Constructor for ruffs_spectra class - a dataframe with spectral data and metadata
#' 
#' @param wavelength Numeric vector of wavelengths (nm)
#' @param values Numeric vector of spectral values (reflectance, radiance, etc.)
#' @param value_name Name for the value column (e.g., "reflectance", "radiance")
#' @param metadata Named list of metadata to attach as attributes
#' @return A ruffs_spectra object
#' @keywords internal
new_ruffs_spectra <- function(wavelength, values, value_name = "reflectance", 
                              metadata = list()) {
  # Create base dataframe
  df <- data.frame(
    wavelength = wavelength,
    value = values,
    stringsAsFactors = FALSE
  )
  names(df)[2] <- value_name
  
  # Attach all metadata as attributes
  for (name in names(metadata)) {
    attr(df, name) <- metadata[[name]]
  }
  
  # Set class
  class(df) <- c("ruffs_spectra", "data.frame")
  
  return(df)
}

#' Create a RUFFS spectra collection
#' 
#' Constructor for ruffs_spectra_collection class - a list of ruffs_spectra objects
#' 
#' @param spectra_list List of ruffs_spectra objects
#' @return A ruffs_spectra_collection object
#' @keywords internal
new_ruffs_spectra_collection <- function(spectra_list) {
  # Validate all elements are ruffs_spectra objects
  is_valid <- sapply(spectra_list, inherits, "ruffs_spectra")
  if (!all(is_valid)) {
    stop("All elements must be ruffs_spectra objects")
  }
  
  class(spectra_list) <- c("ruffs_spectra_collection", "list")
  return(spectra_list)
}

#' Print method for ruffs_spectra
#' @export
print.ruffs_spectra <- function(x, ...) {
  cat("RUFFS Spectral Data\n")
  cat("===================\n")
  
  # Get data type and column name
  dtype <- attr(x, "data_type")
  if (is.null(dtype)) dtype <- "unknown"
  data_col <- names(x)[names(x) != "wavelength"][1]
  
  cat(sprintf("Data type: %s\n", dtype))
  cat(sprintf("Wavelength range: %.2f - %.2f nm\n", 
              min(x$wavelength, na.rm = TRUE), 
              max(x$wavelength, na.rm = TRUE)))
  cat(sprintf("Number of bands: %d\n", nrow(x)))
  
  if (nrow(x) > 1) {
    cat(sprintf("Spectral resolution: %.2f nm (mean)\n", 
                mean(diff(x$wavelength), na.rm = TRUE)))
  }
  
  # Print key metadata
  if (!is.null(attr(x, "instrument_type"))) {
    cat(sprintf("Instrument: %s", attr(x, "instrument_type")))
    if (!is.null(attr(x, "instrument_model"))) {
      cat(sprintf(" (%s)", attr(x, "instrument_model")))
    }
    cat("\n")
  }
  
  if (!is.null(attr(x, "filename"))) {
    cat(sprintf("File: %s\n", attr(x, "filename")))
  }
  
  if (!is.null(attr(x, "acquisition_time"))) {
    cat(sprintf("Acquired: %s\n", attr(x, "acquisition_time")))
  }
  
  cat("\nSpectral data:\n")
  print(head(x, 5))
  if (nrow(x) > 5) {
    cat(sprintf("... with %d more rows\n", nrow(x) - 5))
  }
  
  invisible(x)
}

#' Print method for ruffs_spectra_collection
#' @export
print.ruffs_spectra_collection <- function(x, ...) {
  cat("RUFFS Spectral Data Collection\n")
  cat("===============================\n")
  cat(sprintf("Number of spectra: %d\n\n", length(x)))
  
  if (length(x) > 0) {
    cat("Summary of spectra:\n")
    n_display <- min(10, length(x))
    
    for (i in 1:n_display) {
      name <- names(x)[i]
      if (is.null(name) || name == "") name <- sprintf("[[%d]]", i)
      
      spec <- x[[i]]
      cat(sprintf("  %s: %d bands (%.1f-%.1f nm)",
                  name,
                  nrow(spec),
                  min(spec$wavelength, na.rm = TRUE),
                  max(spec$wavelength, na.rm = TRUE)))
      
      if (!is.null(attr(spec, "instrument_type"))) {
        cat(sprintf(" [%s]", attr(spec, "instrument_type")))
      }
      cat("\n")
    }
    
    if (length(x) > 10) {
      cat(sprintf("  ... and %d more\n", length(x) - 10))
    }
  }
  
  invisible(x)
}

#' Subset method for ruffs_spectra that preserves attributes
#' @export
`[.ruffs_spectra` <- function(x, i, j, drop = FALSE) {
  result <- NextMethod()
  
  if (is.data.frame(result)) {
    # Preserve RUFFS-specific attributes
    attrs_to_keep <- attributes(x)
    attrs_to_keep <- attrs_to_keep[!names(attrs_to_keep) %in% 
                                     c("names", "row.names", "class")]
    attributes(result) <- c(attributes(result), attrs_to_keep)
    class(result) <- class(x)
  }
  
  return(result)
}

#' Extract metadata from RUFFS spectra
#' 
#' @param spectra A ruffs_spectra object
#' @param field Optional specific metadata field to extract
#' @return List of all metadata, or specific field value
#' @export
#' @examples
#' \dontrun{
#' spec <- asdreader("test.asd")
#' get_metadata(spec)
#' get_metadata(spec, "instrument_type")
#' }
get_metadata <- function(spectra, field = NULL) {
  if (!inherits(spectra, "ruffs_spectra")) {
    stop("Input must be a ruffs_spectra object")
  }
  
  attrs <- attributes(spectra)
  attrs <- attrs[!names(attrs) %in% c("names", "row.names", "class")]
  
  if (!is.null(field)) {
    if (field %in% names(attrs)) {
      return(attrs[[field]])
    } else {
      warning(sprintf("Field '%s' not found in metadata", field))
      return(NULL)
    }
  }
  
  return(attrs)
}

#' Apply function to collection of spectra
#' 
#' @param spectra_collection A ruffs_spectra_collection object
#' @param fun Function to apply to each spectrum
#' @param ... Additional arguments passed to fun
#' @return A ruffs_spectra_collection with results
#' @export
#' @examples
#' \dontrun{
#' interpolated <- apply_to_collection(specs, interpolation, 
#'                                      target_wavelengths = seq(350, 2400, 1))
#' }
apply_to_collection <- function(spectra_collection, fun, ...) {
  if (!inherits(spectra_collection, "ruffs_spectra_collection")) {
    stop("Input must be a ruffs_spectra_collection object")
  }
  
  result <- lapply(spectra_collection, fun, ...)
  class(result) <- class(spectra_collection)
  
  return(result)
}