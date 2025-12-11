#' Convert RUFFS Collection to Tidyverse Tibble
#'
#' Converts a RUFFS collection into a long-format tibble suitable for use with 
#' tidyverse packages like dplyr and ggplot2. Supports both ruffs_spectra_summary 
#' and ruffs_spectra_collection objects.
#'
#' @param collection A ruffs_spectra_summary or ruffs_spectra_collection object
#' @param id_column Character string specifying the name for the identifier 
#'   column (default: "spectrum_id")
#'
#' @return A tibble with columns for the identifier and spectral data. For 
#'   ruffs_spectra_summary: wavelength, mean, sd, se, n. For 
#'   ruffs_spectra_collection: wavelength, reflectance (or other data columns).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert summary collection to tibble
#' summary_tibble <- to_tibble(my_summary_collection)
#' 
#' # Convert spectral collection to tibble
#' spectra_tibble <- to_tibble(my_spectral_collection)
#' 
#' # Use with tidyverse tools
#' library(dplyr)
#' summary_tibble %>%
#'   filter(wavelength >= 400, wavelength <= 700)
#' }
to_tibble <- function(collection, id_column = "spectrum_id") {
  
  # Get all classes
  obj_classes <- class(collection)
  
  # Check object class and dispatch to appropriate method
  # Check for spectra_collection first (more common use case)
  if ("ruffs_spectra_collection" %in% obj_classes) {
    return(.to_tibble_collection(collection, id_column))
  } else if ("ruffs_spectra_summary" %in% obj_classes) {
    return(.to_tibble_summary(collection, id_column))
  } else {
    stop("Input must be a ruffs_spectra_summary or ruffs_spectra_collection object.\n",
         "Found classes: ", paste(obj_classes, collapse = ", "))
  }
}

#' Convert Summary Collection to Tibble (Internal)
#'
#' @param collection A ruffs_spectra_summary collection object
#' @param id_column Character string for identifier column name
#' @return A tibble
#' @noRd
.to_tibble_summary <- function(collection, id_column) {
  # Extract each dataframe and add its name as an ID column
  result <- purrr::map2_dfr(
    collection, 
    names(collection),
    ~ dplyr::mutate(.x, !!id_column := .y)
  ) %>%
    tibble::as_tibble()
  
  return(result)
}

#' Convert Spectral Collection to Tibble (Internal)
#'
#' @param collection A ruffs_spectra_collection object
#' @param id_column Character string for identifier column name
#' @return A tibble
#' @noRd
.to_tibble_collection <- function(collection, id_column) {
  # Extract each dataframe and add its name as an ID column
  result <- purrr::map2_dfr(
    collection, 
    names(collection),
    ~ dplyr::mutate(.x, !!id_column := .y)
  ) %>%
    tibble::as_tibble()
  
  return(result)
}
