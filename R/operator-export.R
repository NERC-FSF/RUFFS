#' Convert RUFFS Collection to Tidyverse Tibble
#'
#' Converts a ruffs_spectra_summary collection into a long-format tibble
#' suitable for use with tidyverse packages like dplyr and ggplot2.
#'
#' @param collection A ruffs_spectra_summary collection object
#' @param id_column Character string specifying the name for the identifier 
#'   column (default: "spectrum_id")
#'
#' @return A tibble with columns for the identifier and all spectral data 
#'   columns in ruffs_spectra_summary collection (wavelength, mean, sd, se, n)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert summary collection to tibble
#' summary_tibble <- to_tibble(my_summary_collection)
#' 
#' # Use with tidyverse tools
#' library(dplyr)
#' summary_tibble %>%
#'   filter(wavelength >= 400, wavelength <= 700)
#' }
to_tibble <- function(collection, id_column = "spectrum_id") {
  if (!inherits(collection, "ruffs_spectra_summary")) {
    stop("Input must be a ruffs_spectra_summary object")
  }
  
  # Extract each dataframe and add its name as an ID column
  result <- purrr::map2_dfr(
    collection, 
    names(collection),
    ~ dplyr::mutate(.x, !!id_column := .y)
  ) %>%
    tibble::as_tibble()
  
  return(result)
}