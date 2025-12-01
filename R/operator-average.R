#' Average Spectra by Group
#'
#' @description
#' Computes the mean spectrum for each group in a grouped spectral collection.
#' This reduces multiple replicate measurements to a single representative
#' spectrum per group.
#'
#' @param spectra A `ruffs_spectra_collection_grouped` object created by
#'   [group_spectra()].
#' @param na.rm Logical. If `TRUE` (default), NA values are removed before
#'   computing means. If `FALSE`, any NA values will propagate to the result.
#'
#' @return A `ruffs_spectra_collection` object containing one averaged spectrum
#'   per group. Each spectrum is a dataframe with columns:
#'   \itemize{
#'     \item `wavelength`: Wavelength values (nm)
#'     \item `spectra`: Mean spectral values across all group members
#'   }
#'   The collection retains metadata from the input (instrument, data_type,
#'   wavelength_range) and adds:
#'   \itemize{
#'     \item `averaged`: Logical flag set to `TRUE`
#'     \item `group_sizes`: Named integer vector with number of spectra
#'       averaged in each group
#'   }
#'
#' @details
#' The function assumes all spectra within a group share identical wavelength
#' grids (as would be the case after using [interpolate_spectra()]).
#'
#' Group names become the names of the averaged spectra in the output
#' collection. For example, if groups are "site1_oak" and "site1_pine", the
#' output collection will contain spectra named "site1_oak" and "site1_pine".
#'
#' @examples
#' \dontrun{
#' # Process and group spectra
#' spectra <- svcreader("path/to/files") |>
#'   interpolate_spectra() |>
#'   stitch_spectra() |>
#'   correct_spectra() |>
#'   group_spectra(separator = "_", position = 2)
#'
#' # Compute group averages
#' averaged <- average_spectra(spectra)
#'
#' # View result
#' print(averaged)
#' names(averaged$spectra)  # Group names
#'
#' # Access individual averaged spectrum
#' averaged$spectra[["site1_oak"]]
#'
#' # Check how many spectra were averaged per group
#' attr(averaged, "group_sizes")
#' }
#'
#' @seealso [group_spectra()] for creating grouped collections,
#'   [summarise_spectra()] for computing additional statistics (SD, SE, etc.)
#'
#' @export
average_spectra <- function(spectra, na.rm = TRUE) {


  # Input validation
  .validate_average_inputs(spectra)

  # Get group information
  group_names <- get_group_names(spectra)
  group_assignments <- get_group_assignments(spectra)

  # Compute average for each group
  averaged_list <- lapply(group_names, function(grp) {
    .compute_group_mean(spectra, grp, na.rm = na.rm)
  })
  names(averaged_list) <- group_names

  # Track group sizes
  group_sizes <- .get_group_sizes(spectra, group_names)

  # Create output collection
  result <- .create_averaged_collection(
    averaged_list = averaged_list,
    group_sizes = group_sizes,
    original_spectra = spectra

  )

  # Report results
  .report_averaging(group_names, group_sizes)

  return(result)
}


#' Validate inputs for average_spectra
#' @noRd
.validate_average_inputs <- function(spectra) {

  if (!inherits(spectra, "ruffs_spectra_collection_grouped")) {
    stop(
      "Input must be a grouped spectra collection.\n",
      "Use group_spectra() first to group your collection.\n",
      "Received class: ", class(spectra)[1],
      call. = FALSE
    )
  }

  # Check that groups exist
  group_names <- attr(spectra, "group_names")
  if (is.null(group_names) || length(group_names) == 0) {
    stop("No groups found in the collection.", call. = FALSE)
  }
}


#' Compute mean spectrum for a single group
#' @noRd
.compute_group_mean <- function(spectra, group_name, na.rm = TRUE) {

  # Extract spectra for this group
  group_data <- get_spectra_by_group(spectra, group_name)
  spectra_list <- group_data$spectra

  # Get first spectrum to determine structure
  first_spectrum <- spectra_list[[1]]

  # Determine the spectral value column name
  # Common names: "spectra", "reflectance", "radiance", "value", "dn"
  spec_col <- .get_spectral_column(first_spectrum)

  # Get wavelengths
  wavelengths <- first_spectrum$wavelength

  # Handle single spectrum case
  if (length(spectra_list) == 1) {
    return(data.frame(
      wavelength = wavelengths,
      spectra = first_spectrum[[spec_col]]
    ))
  }

  # Extract spectral values into matrix (wavelengths × spectra)
  # Use do.call + cbind to ensure matrix output
  values_list <- lapply(spectra_list, function(s) s[[spec_col]])
  values_matrix <- do.call(cbind, values_list)

  # Compute row means
  mean_values <- rowMeans(values_matrix, na.rm = na.rm)

  data.frame(
    wavelength = wavelengths,
    spectra = mean_values
  )
}


#' Identify the spectral value column in a spectrum dataframe
#' @noRd
.get_spectral_column <- function(spectrum_df) {

  col_names <- names(spectrum_df)


  # Look for common spectral column names (excluding wavelength)
  candidates <- c("spectra", "reflectance", "radiance", "value", "dn", "raw_dn")

  for (candidate in candidates) {
    if (candidate %in% col_names) {
      return(candidate)
    }
  }

  # Fallback: use second column (assuming first is wavelength)
  non_wl_cols <- setdiff(col_names, c("wavelength", "Wavelength", "wl", "WL"))

  if (length(non_wl_cols) >= 1) {
    return(non_wl_cols[1])
  }

  stop(
    "Cannot identify spectral value column.\n",
    "Available columns: ", paste(col_names, collapse = ", "),
    call. = FALSE
  )
}


#' Get the number of spectra in each group
#' @noRd
.get_group_sizes <- function(spectra, group_names) {

  group_assignments <- get_group_assignments(spectra)

  sizes <- vapply(group_names, function(grp) {
    sum(group_assignments == grp)
  }, integer(1))

  return(sizes)
}


#' Create the output averaged collection
#' @noRd
.create_averaged_collection <- function(averaged_list, group_sizes,
                                        original_spectra) {

  # Build collection structure
  result <- list(spectra = averaged_list)

  # Set class
  class(result) <- "ruffs_spectra_collection"

  # Preserve original metadata
  attr(result, "instrument") <- attr(original_spectra, "instrument")
  attr(result, "data_type") <- attr(original_spectra, "data_type")
  attr(result, "wavelength_range") <- attr(original_spectra, "wavelength_range")

  # Add averaging metadata
  attr(result, "averaged") <- TRUE
  attr(result, "group_sizes") <- group_sizes

  return(result)
}


#' Report averaging results to console
#' @noRd
.report_averaging <- function(group_names, group_sizes) {

  n_groups <- length(group_names)
  total_spectra <- sum(group_sizes)

  message(
    "Averaged ", total_spectra, " spectra into ", n_groups, " group mean(s):"
  )

  for (i in seq_along(group_names)) {
    message("  - '", group_names[i], "': n = ", group_sizes[i])
  }
}
