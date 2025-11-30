#' Interpolate spectral data to target wavelengths
#' 
#' Interpolates spectral data to a new set of wavelengths.
#' Handles irregular spacing and non-monotonic sequences in the original data.
#' 
#' @param spectra A ruffs_spectra object
#' @param target_wavelengths Numeric vector of target wavelengths (in nm)
#' @param method Interpolation method. Options:
#'   \describe{
#'     \item{"linear"}{Linear interpolation (default)}
#'     \item{"spline"}{Cubic spline interpolation}
#'     \item{"constant"}{Nearest-neighbor (constant) interpolation}
#'   }
#' @return A ruffs_spectra object with interpolated data
#' @export
#' @examples
#' \dontrun{
#' # Read SVC data
#' spec <- svcreader("sample.sig")
#' 
#' # Interpolate to 1nm spacing
#' spec_interp <- interpolation(spec, seq(400, 2400, 1))
#' 
#' # Use spline interpolation
#' spec_spline <- interpolation(spec, seq(400, 2400, 1), method = "spline")
#' }
interpolation <- function(spectra, target_wavelengths, method = "linear") {
  
  # Input validation
  if (!inherits(spectra, "ruffs_spectra")) {
    stop("spectra must be a ruffs_spectra object")
  }
  
  if (!is.numeric(target_wavelengths)) {
    stop("target_wavelengths must be a numeric vector")
  }
  
  if (length(target_wavelengths) == 0) {
    stop("target_wavelengths must contain at least one value")
  }
  
  # Validate method
  valid_methods <- c("linear", "spline", "constant")
  if (!method %in% valid_methods) {
    stop(sprintf("method must be one of: %s", paste(valid_methods, collapse = ", ")))
  }
  
  # Check for wavelength column
  if (!"wavelength" %in% names(spectra)) {
    stop("spectra must contain a 'wavelength' column")
  }
  
  # Get the data column (not wavelength)
  data_cols <- setdiff(names(spectra), "wavelength")
  if (length(data_cols) == 0) {
    stop("spectra must contain at least one data column in addition to wavelength")
  }
  
  # Use the first data column
  data_col <- data_cols[1]
  
  # Check for non-monotonic wavelengths
  if (is.unsorted(spectra$wavelength)) {
    warning("Wavelengths are not monotonically increasing. Sorting data before interpolation.")
    sort_idx <- order(spectra$wavelength)
    spectra <- spectra[sort_idx, ]
  }
  
  # Remove any duplicate wavelengths (keep first occurrence)
  if (any(duplicated(spectra$wavelength))) {
    warning("Duplicate wavelengths detected. Keeping first occurrence of each wavelength.")
    spectra <- spectra[!duplicated(spectra$wavelength), ]
  }
  
  # Perform interpolation based on method
  if (method == "linear") {
    interp_values <- approx(
      x = spectra$wavelength,
      y = spectra[[data_col]],
      xout = target_wavelengths,
      method = "linear",
      rule = 2  # Use nearest value for extrapolation
    )$y
    
  } else if (method == "spline") {
    # Use spline interpolation
    interp_values <- spline(
      x = spectra$wavelength,
      y = spectra[[data_col]],
      xout = target_wavelengths,
      method = "natural"
    )$y
    
  } else if (method == "constant") {
    # Nearest neighbor interpolation
    interp_values <- approx(
      x = spectra$wavelength,
      y = spectra[[data_col]],
      xout = target_wavelengths,
      method = "constant",
      rule = 2
    )$y
  }
  
  # Get current metadata
  current_metadata <- get_metadata(spectra)
  
  # Create output with same metadata
  result <- new_ruffs_spectra(
    wavelength = target_wavelengths,
    values = interp_values,
    value_name = data_col,
    metadata = current_metadata
  )
  
  # Add processing note to metadata
  processing_note <- sprintf(
    "Interpolated from %d to %d wavelengths using %s method (range: %.2f-%.2f nm)",
    nrow(spectra),
    length(target_wavelengths),
    method,
    min(target_wavelengths),
    max(target_wavelengths)
  )
  
  attr(result, "processing_history") <- c(
    attr(spectra, "processing_history"),
    processing_note
  )
  
  return(result)
}