#' Jump Correction for Spectral Data
#'
#' Corrects for discontinuities (jumps) at detector boundaries after overlap
#' stitching. After overlapping wavelength regions are removed by 
#' \code{overlap_stitching()}, small step artifacts can remain at the detector
#' junctions. This function applies additive corrections to create smooth
#' transitions between detector regions.
#'
#' @param data A data frame containing spectral data with a 'wavelength' column
#'   and one or more measurement columns, or a ruffs_spectra object.
#' @param splices Numeric vector of wavelength values indicating detector 
#'   boundaries where jumps may occur. For ASD instruments, typically 
#'   c(1000, 1800). For SVC instruments, c(990, 1850).
#' @param reference Integer indicating which detector segment to use as the
#'   reference (1-based indexing). This segment remains unchanged while others
#'   are adjusted to match it. Default is 2 (middle detector for 3-detector systems).
#' @param method Character string specifying the correction method. Currently
#'   only "additive" is supported.
#' @param measurement_cols Character vector of column names to apply correction to.
#'   If NULL (default), applies to all numeric columns except 'wavelength'.
#'
#' @return A data frame or ruffs_spectra object with jump-corrected measurements.
#'   The structure and class of the input are preserved.
#'
#' @details
#' The additive correction method works by:
#' \enumerate{
#'   \item Dividing the spectrum into segments based on splice points
#'   \item Keeping the reference segment unchanged
#'   \item For segments to the right: shifting each so its first value matches
#'         the last value of the previous segment
#'   \item For segments to the left: shifting each so its last value matches
#'         the first value of the next segment
#' }
#'
#' This creates a continuous spectrum by removing step discontinuities while
#' preserving the relative spectral shape within each detector region.
#'
#' **Important:** This function should be applied AFTER \code{overlap_stitching()}
#' has removed overlapping wavelength regions. If overlaps still exist, the
#' correction may not work as intended.
#'
#' @section Choosing Splice Points:
#' Splice points should match the detector boundaries of your instrument:
#' \itemize{
#'   \item ASD instruments: typically 1000 nm (VNIR/SWIR1) and 1800 nm (SWIR1/SWIR2)
#'   \item SVC instruments: typically 990 nm and 1800 nm
#'   \item PSR instruments: typically 1000 nm and 1800 nm
#' }
#'
#' @section Choosing Reference Segment:
#' The reference segment is the detector region that remains unchanged:
#' \itemize{
#'   \item For 3-detector systems: reference = 2 (middle detector, usually SWIR1)
#'   \item For 2-detector systems: reference = 1 or 2 depending on which is more reliable
#' }
#'
#' @examples
#' \dontrun{
#' # Single spectrum processing
#' data <- svcreader("field_data.sig", data_type = "reflectance")
#' data_interp <- interpolation(data, seq(350, 2500, 1))
#' data_stitched <- overlap_stitching(data_interp, method = "average")
#' data_corrected <- jump_correct(data_stitched, splices = c(990, 1850), reference = 2)
#'
#' # Collection processing (RUFFS pattern)
#' svcdata <- svcreader(directory = "data/", data_type = "reflectance", recursive = TRUE)
#' specs_interpolation <- apply_to_collection(svcdata, interpolation, seq(350, 2500, 1))
#' specs_stitched <- apply_to_collection(specs_interpolation, overlap_stitching, method = "average")
#' specs_corrected <- apply_to_collection(specs_stitched, jump_correct, 
#'                                        splices = c(990, 1850), reference = 2)
#'
#' # For ASD data
#' asd_data <- asdreader("spectrum.asd", data_type = "reflectance")
#' asd_corrected <- asd_data %>%
#'   interpolation(seq(350, 2500, 1)) %>%
#'   overlap_stitching(method = "average") %>%
#'   jump_correct(splices = c(1000, 1800), reference = 2)
#'
#' # Apply to specific columns only
#' data_corrected <- jump_correct(data_stitched, 
#'                                 splices = c(990, 1850),
#'                                 reference = 2,
#'                                 measurement_cols = c("reflectance", "radiance"))
#' }
#'
#' @export
jump_correct <- function(data, splices, reference = 2, method = "additive", 
                         measurement_cols = NULL) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame or ruffs_spectra object")
  }
  
  if (!"wavelength" %in% names(data)) {
    stop("Data must contain a 'wavelength' column")
  }
  
  if (!is.numeric(splices) || length(splices) == 0) {
    stop("splices must be a numeric vector with at least one value")
  }
  
  if (!is.numeric(reference) || length(reference) != 1 || reference < 1) {
    stop("reference must be a positive integer")
  }
  
  if (reference > (length(splices) + 1)) {
    stop("reference must be <= number of segments (", length(splices) + 1, ")")
  }
  
  if (method != "additive") {
    stop("Currently only 'additive' method is supported")
  }
  
  # Check for non-monotonic wavelengths
  if (any(diff(data$wavelength) <= 0)) {
    stop("Wavelength column must be monotonically increasing. ",
         "This suggests overlapping regions still exist. ",
         "Apply overlap_stitching() before jump_correct().")
  }
  
  # Store original attributes
  original_attrs <- attributes(data)
  original_class <- class(data)
  
  # Determine measurement columns
  if (is.null(measurement_cols)) {
    measurement_cols <- setdiff(names(data), "wavelength")
    measurement_cols <- measurement_cols[sapply(data[measurement_cols], is.numeric)]
  } else {
    # Validate specified columns
    missing_cols <- setdiff(measurement_cols, names(data))
    if (length(missing_cols) > 0) {
      stop("Specified measurement columns not found: ", 
           paste(missing_cols, collapse = ", "))
    }
    non_numeric <- measurement_cols[!sapply(data[measurement_cols], is.numeric)]
    if (length(non_numeric) > 0) {
      stop("Specified measurement columns must be numeric: ",
           paste(non_numeric, collapse = ", "))
    }
  }
  
  if (length(measurement_cols) == 0) {
    warning("No numeric measurement columns found. Returning data unchanged.")
    return(data)
  }
  
  # Apply correction to each measurement column
  for (col in measurement_cols) {
    data[[col]] <- .jump_correct_additive(
      wavelength = data$wavelength,
      measurement = data[[col]],
      splices = splices,
      reference = reference
    )
  }
  
  # Restore class and attributes
  class(data) <- original_class
  for (attr_name in setdiff(names(original_attrs), c("names", "row.names", "class"))) {
    attr(data, attr_name) <- original_attrs[[attr_name]]
  }
  
  return(data)
}


#' Internal function for additive jump correction
#'
#' @param wavelength Numeric vector of wavelengths
#' @param measurement Numeric vector of measurements
#' @param splices Numeric vector of splice points
#' @param reference Integer indicating reference segment (1-based)
#'
#' @return Numeric vector of corrected measurements
#' @keywords internal
.jump_correct_additive <- function(wavelength, measurement, splices, reference) {
  
  # Create a copy to modify
  corrected <- measurement
  
  # Assign each wavelength to a segment
  segment_id <- .get_segment_id(wavelength, splices)
  
  # Get unique segments in order
  unique_segments <- sort(unique(segment_id))
  n_segments <- length(unique_segments)
  
  # Validate reference
  if (reference > n_segments) {
    stop("Reference segment (", reference, ") exceeds number of segments (", 
         n_segments, ")")
  }
  
  # Helper function to apply translation
  .translate_segment <- function(ref_idx, mov_idx, to_right = TRUE) {
    ref_segment <- segment_id == unique_segments[ref_idx]
    mov_segment <- segment_id == unique_segments[mov_idx]
    
    ref_values <- corrected[ref_segment]
    mov_values <- corrected[mov_segment]
    
    # Calculate offset
    if (to_right) {
      # Match first value of moving segment to last value of reference
      offset <- ref_values[length(ref_values)] - mov_values[1]
    } else {
      # Match last value of moving segment to first value of reference
      offset <- ref_values[1] - mov_values[length(mov_values)]
    }
    
    # Apply offset to moving segment
    corrected[mov_segment] <<- mov_values + offset
  }
  
  # Process segments to the right of reference
  if (reference < n_segments) {
    for (i in reference:(n_segments - 1)) {
      .translate_segment(i, i + 1, to_right = TRUE)
    }
  }
  
  # Process segments to the left of reference
  if (reference > 1) {
    for (i in reference:2) {
      .translate_segment(i, i - 1, to_right = FALSE)
    }
  }
  
  return(corrected)
}


#' Assign wavelengths to segments based on splice points
#'
#' @param wavelength Numeric vector of wavelengths
#' @param splices Numeric vector of splice points
#'
#' @return Integer vector of segment IDs (1-based)
#' @keywords internal
.get_segment_id <- function(wavelength, splices) {
  n <- length(wavelength)
  segment_id <- integer(n)
  
  for (i in seq_along(wavelength)) {
    wl <- wavelength[i]
    segment <- 1
    
    for (splice in splices) {
      if (wl > splice) {
        segment <- segment + 1
      } else {
        break
      }
    }
    
    segment_id[i] <- segment
  }
  
  return(segment_id)
}