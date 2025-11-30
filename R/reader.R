#' Universal spectrometer file reader
#' 
#' Automatically detects instrument type and reads spectral data from ASD, SVC, or PSR files.
#' This is the main entry point for reading spectral data in RUFFS.
#' 
#' @param file_path Character vector of file path(s) to spectral data file(s).
#'   If NULL, reads all compatible files from current working directory.
#' @param instrument_type Optional instrument type specification. If NULL (default),
#'   will be auto-detected from file extension. Options: "asd", "svc", "psr"
#' @param data_type Type of data to extract. Default is "reflectance".
#'   See individual reader functions for instrument-specific options.
#' @return A ruffs_spectra object (single file) or ruffs_spectra_collection (multiple files)
#' @export
#' @examples
#' \dontrun{
#' # Read single file with auto-detection
#' spec <- reader("sample.asd")
#' 
#' # Read all spectral files in current directory
#' specs <- reader()
#' 
#' # Read multiple files
#' specs <- reader(c("file1.sig", "file2.sig", "file3.sig"))
#' 
#' # Force specific reader and data type
#' spec <- reader("data.txt", instrument_type = "svc", data_type = "radiance_target")
#' 
#' # Read radiance from ASD files
#' rad <- reader("sample.asd", data_type = "radiance")
#' }
reader <- function(file_path = NULL, instrument_type = NULL, data_type = "reflectance") {
  
  # Handle NULL file_path - try to auto-detect files in current directory
  if (is.null(file_path)) {
    if (is.null(instrument_type)) {
      # Look for all supported file types
      asd_files <- list.files(pattern = "\\.asd$", ignore.case = TRUE, full.names = TRUE)
      sig_files <- list.files(pattern = "\\.sig$", ignore.case = TRUE, full.names = TRUE)
      sed_files <- list.files(pattern = "\\.sed$", ignore.case = TRUE, full.names = TRUE)
      
      all_files <- c(asd_files, sig_files, sed_files)
      
      if (length(all_files) == 0) {
        stop("No supported spectral files (.asd, .sig, .sed) found in current directory")
      }
      
      # Check if multiple types exist
      types_found <- c(
        if(length(asd_files) > 0) "ASD" else NULL,
        if(length(sig_files) > 0) "SVC" else NULL,
        if(length(sed_files) > 0) "PSR" else NULL
      )
      
      if (length(types_found) > 1) {
        stop(sprintf("Multiple instrument types found (%s). Please specify instrument_type.",
                     paste(types_found, collapse = ", ")))
      }
      
      file_path <- all_files
      message(sprintf("Found %d %s file(s) in current directory", 
                      length(file_path), types_found[1]))
    } else {
      # instrument_type specified but no file_path - let individual reader handle it
      instrument_type <- tolower(instrument_type)
      
      if (!instrument_type %in% c("asd", "svc", "psr")) {
        stop(sprintf("Invalid instrument_type '%s'. Must be one of: asd, svc, psr",
                     instrument_type))
      }
      
      result <- switch(instrument_type,
                       "asd" = asdreader(NULL, data_type = data_type),
                       "svc" = svcreader(NULL, data_type = data_type),
                       "psr" = psrreader(NULL, data_type = data_type)
      )
      
      return(result)
    }
  }
  
  # Validate file_path
  if (!is.character(file_path) || length(file_path) == 0) {
    stop("file_path must be a non-empty character vector or NULL")
  }
  
  # Check all files exist
  missing_files <- file_path[!file.exists(file_path)]
  if (length(missing_files) > 0) {
    stop(sprintf("File(s) not found: %s", 
                 paste(missing_files, collapse = ", ")))
  }
  
  # Auto-detect instrument type if not specified
  if (is.null(instrument_type)) {
    extensions <- tolower(tools::file_ext(file_path))
    
    # Check if all files have same extension
    if (length(unique(extensions)) > 1) {
      stop("Multiple file types detected. Please specify instrument_type explicitly.")
    }
    
    extension <- extensions[1]
    
    instrument_type <- switch(extension,
                              "asd" = "asd",
                              "sig" = "svc",
                              "sed" = "psr",
                              stop(sprintf("Unrecognized file extension: .%s. Please specify instrument_type.", 
                                           extension))
    )
    
    message(sprintf("Auto-detected instrument type: %s", toupper(instrument_type)))
  }
  
  # Normalize and validate instrument type
  instrument_type <- tolower(instrument_type)
  
  if (!instrument_type %in% c("asd", "svc", "psr")) {
    stop(sprintf("Invalid instrument_type '%s'. Must be one of: asd, svc, psr",
                 instrument_type))
  }
  
  # Call appropriate reader with data_type
  result <- switch(instrument_type,
                   "asd" = asdreader(file_path, data_type = data_type),
                   "svc" = svcreader(file_path, data_type = data_type),
                   "psr" = psrreader(file_path, data_type = data_type),
                   stop("Unexpected error in reader dispatch")
  )
  
  return(result)
}