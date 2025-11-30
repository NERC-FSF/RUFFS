#' Stitch overlapping detector regions
#' 
#' Resolves overlapping wavelength regions between detectors where overlap removal has not occured. 
#' 
#' @param spectra A ruffs_spectra object (typically from SVC)
#' @param method Stitching method for overlap regions. Options:
#'   \describe{
#'     \item{"average"}{Average values in overlap region (default)}
#'     \item{"first"}{Use values from first detector}
#'     \item{"second"}{Use values from second detector}
#'     \item{"linear_blend"}{Linear interpolation between detectors}
#'   }
#' @param overlap1 Numeric vector of length 2 defining VNIR-SWIR1 overlap region (nm).
#'   If NULL, uses metadata from spectra or defaults to c(990, 1010)
#' @param overlap2 Numeric vector of length 2 defining SWIR1-SWIR2 overlap region (nm).
#'   If NULL, uses metadata from spectra or defaults to c(1790, 1810)
#' @return A ruffs_spectra object with stitched data (no duplicate wavelengths)
#' @export
#' @examples
#' \dontrun{
#' # Read SVC data
#' spec <- svcreader("sample.sig")
#' 
#' # Stitch overlaps using averaging
#' spec_stitched <- overlap_stitching(spec)
#' 
#' # Use first detector in overlaps
#' spec_first <- overlap_stitching(spec, method = "first")
#' 
#' # Custom overlap regions
#' spec_custom <- overlap_stitching(spec, 
#'                                   overlap1 = c(995, 1005),
#'                                   overlap2 = c(1795, 1805))
#' }
overlap_stitching <- function(spectra, method = "average", 
                              overlap1 = NULL, overlap2 = NULL) {
  
  # Input validation
  if (!inherits(spectra, "ruffs_spectra")) {
    stop("spectra must be a ruffs_spectra object")
  }
  
  # Validate method
  valid_methods <- c("average", "first", "second", "linear_blend")
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
    stop("spectra must contain at least one data column")
  }
  data_col <- data_cols[1]
  
  # Get overlap regions from metadata or use provided/default values
  if (is.null(overlap1)) {
    overlap1 <- attr(spectra, "overlap_vnir_swir1")
    if (is.null(overlap1)) {
      overlap1 <- c(990, 1010)
      message("No overlap metadata found. Using default VNIR-SWIR1 overlap: 990-1010 nm")
    }
  }
  
  if (is.null(overlap2)) {
    overlap2 <- attr(spectra, "overlap_swir1_swir2")
    if (is.null(overlap2)) {
      overlap2 <- c(1790, 1810)
      message("No overlap metadata found. Using default SWIR1-SWIR2 overlap: 1790-1810 nm")
    }
  }
  
  # Identify overlap regions
  overlap1_idx <- which(spectra$wavelength >= overlap1[1] & 
                          spectra$wavelength <= overlap1[2])
  overlap2_idx <- which(spectra$wavelength >= overlap2[1] & 
                          spectra$wavelength <= overlap2[2])
  
  # If no overlaps detected, return original
  if (length(overlap1_idx) == 0 && length(overlap2_idx) == 0) {
    message("No overlap regions detected in the wavelength range. Returning original data.")
    return(spectra)
  }
  
  # Initialize output vectors
  wl <- spectra$wavelength
  values <- spectra[[data_col]]
  
  # Process first overlap region (VNIR-SWIR1)
  if (length(overlap1_idx) > 0) {
    values <- .stitch_overlap_region(
      wavelengths = wl,
      values = values,
      overlap_idx = overlap1_idx,
      overlap_range = overlap1,
      method = method,
      region_name = "VNIR-SWIR1"
    )
  }
  
  # Process second overlap region (SWIR1-SWIR2)
  if (length(overlap2_idx) > 0) {
    values <- .stitch_overlap_region(
      wavelengths = wl,
      values = values,
      overlap_idx = overlap2_idx,
      overlap_range = overlap2,
      method = method,
      region_name = "SWIR1-SWIR2"
    )
  }
  
  # Remove duplicate wavelengths (keep stitched values)
  unique_idx <- !duplicated(wl)
  wl_stitched <- wl[unique_idx]
  values_stitched <- values[unique_idx]
  
  # Get current metadata
  current_metadata <- get_metadata(spectra)
  
  # Create output
  result <- new_ruffs_spectra(
    wavelength = wl_stitched,
    values = values_stitched,
    value_name = data_col,
    metadata = current_metadata
  )
  
  # Add processing note
  processing_note <- sprintf(
    "Overlap stitching applied using '%s' method. Regions: [%.1f-%.1f nm, %.1f-%.1f nm]. Removed %d duplicate wavelengths.",
    method,
    overlap1[1], overlap1[2],
    overlap2[1], overlap2[2],
    sum(!unique_idx)
  )
  
  attr(result, "processing_history") <- c(
    attr(spectra, "processing_history"),
    processing_note
  )
  
  return(result)
}

# Helper function to stitch a single overlap region
.stitch_overlap_region <- function(wavelengths, values, overlap_idx, 
                                   overlap_range, method, region_name) {
  
  if (length(overlap_idx) == 0) {
    return(values)
  }
  
  # Find the indices just before and after the overlap
  pre_overlap_idx <- max(which(wavelengths < overlap_range[1]))
  post_overlap_idx <- min(which(wavelengths > overlap_range[2]))
  
  # Extract overlap values
  overlap_wl <- wavelengths[overlap_idx]
  overlap_vals <- values[overlap_idx]
  
  # Apply stitching method
  if (method == "average") {
    # Group by wavelength and average
    unique_wl <- unique(overlap_wl)
    stitched_vals <- sapply(unique_wl, function(w) {
      mean(overlap_vals[overlap_wl == w], na.rm = TRUE)
    })
    
    # Replace overlap region with averaged values
    values[overlap_idx[1:length(unique_wl)]] <- stitched_vals
    
  } else if (method == "first") {
    # Keep first occurrence of each wavelength
    keep_idx <- !duplicated(overlap_wl)
    values[overlap_idx] <- overlap_vals[keep_idx]
    
  } else if (method == "second") {
    # Keep last occurrence of each wavelength
    keep_idx <- !duplicated(overlap_wl, fromLast = TRUE)
    values[overlap_idx] <- overlap_vals[keep_idx]
    
  } else if (method == "linear_blend") {
    # Linear interpolation from first to second detector
    overlap_length <- length(overlap_idx)
    if (overlap_length > 1) {
      blend_weights <- seq(0, 1, length.out = overlap_length)
      
      # Get first and last values in overlap
      first_vals <- values[overlap_idx[1]]
      last_vals <- values[overlap_idx[overlap_length]]
      
      # Blend
      for (i in 1:overlap_length) {
        values[overlap_idx[i]] <- (1 - blend_weights[i]) * first_vals + 
          blend_weights[i] * last_vals
      }
    }
  }
  
  return(values)
}