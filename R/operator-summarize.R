#' Summarise Spectra by Group
#'
#' @description
#' Computes summary statistics (mean, standard deviation, standard error,
#' minimum, maximum, and sample size) for each group in a grouped spectral
#' collection.
#'
#' @param spectra A `ruffs_spectra_collection_grouped` object created by
#'   [group_spectra()].
#' @param stats Character vector specifying which statistics to compute.
#'   Options are: `"mean"`, `"sd"`, `"se"`, `"min"`, `"max"`, `"median"`,
#'   `"n"`. Default is `c("mean", "sd", "se", "n")`.
#' @param na.rm Logical. If `TRUE` (default), NA values are removed before
#'   computing statistics.
#'
#' @return A `ruffs_spectra_summary` object, which is a named list where each
#'   element corresponds to a group. Each group element is a dataframe with:
#'   \itemize{
#'     \item `wavelength`: Wavelength values (nm)
#'     \item Additional columns for each requested statistic (mean, sd, etc.)
#'   }
#'   The object also contains attributes:
#'   \itemize{
#'     \item `group_names`: Character vector of group names
#'     \item `stats_computed`: Character vector of statistics that were computed
#'     \item `instrument`, `data_type`, `wavelength_range`: Preserved from input
#'   }
#'
#' @details
#' Statistics available:
#' \describe{
#'   \item{mean}{Arithmetic mean across group members}
#'   \item{sd}{Sample standard deviation}
#'   \item{se}{Standard error of the mean (sd / sqrt(n))}
#'   \item{min}{Minimum value}
#'   \item{max}{Maximum value}
#'   \item{median}{Median value}
#'   \item{n}{Sample size (number of spectra in group)}
#' }
#'
#' The function assumes all spectra within a group share identical wavelength
#' grids (as would be the case after using [interpolate_spectra()]).
#'
#' For groups with only one spectrum, `sd` and `se` will be `NA`.
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
#' # Compute default statistics (mean, sd, se, n)
#' summary <- summarise_spectra(spectra)
#'
#' # Compute all available statistics
#' summary_full <- summarise_spectra(spectra,
#'   stats = c("mean", "sd", "se", "min", "max", "median", "n"))
#'
#' # Access summary for a specific group
#' summary[["site1_oak"]]
#' head(summary[["site1_oak"]])
#' #   wavelength       mean          sd          se n
#' # 1        350 0.05234521 0.002345210 0.001354213 3
#' # 2        351 0.05312453 0.002456321 0.001418324 3
#' # ...
#'
#' # Plot mean with error bars
#' grp_data <- summary[["site1_oak"]]
#' plot(grp_data$wavelength, grp_data$mean, type = "l")
#' # Add ± 1 SE bands
#' lines(grp_data$wavelength, grp_data$mean + grp_data$se, lty = 2)
#' lines(grp_data$wavelength, grp_data$mean - grp_data$se, lty = 2)
#' }
#'
#' @seealso [group_spectra()] for creating grouped collections,
#'   [average_spectra()] for computing just the mean,
#'   [plot_summary()] for visualising summary statistics
#'
#' @export
summarise_spectra <- function(spectra,
                              stats = c("mean", "sd", "se", "n"),
                              na.rm = TRUE) {
  
  # Input validation
  .validate_summary_inputs(spectra, stats)
  
  # Standardise stats names
  stats <- tolower(stats)
  
  # Get group information
  group_names <- get_group_names(spectra)
  
  # Compute summary for each group
  summary_list <- lapply(group_names, function(grp) {
    .compute_group_summary(spectra, grp, stats = stats, na.rm = na.rm)
  })
  names(summary_list) <- group_names
  
  # Create output object
  result <- .create_summary_object(
    summary_list = summary_list,
    group_names = group_names,
    stats = stats,
    original_spectra = spectra
  )
  
  # Report results
  .report_summary(group_names, stats)
  
  return(result)
}


#' Alias for British English spelling
#' @rdname summarise_spectra
#' @export
summarize_spectra <- summarise_spectra


#' Validate inputs for summarise_spectra
#' @noRd
.validate_summary_inputs <- function(spectra, stats) {
  
  # Check spectra class
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
  
  # Validate stats argument
  valid_stats <- c("mean", "sd", "se", "min", "max", "median", "n")
  stats_lower <- tolower(stats)
  invalid <- setdiff(stats_lower, valid_stats)
  
  if (length(invalid) > 0) {
    stop(
      "Invalid statistic(s) requested: ", paste(invalid, collapse = ", "), "\n",
      "Valid options are: ", paste(valid_stats, collapse = ", "),
      call. = FALSE
    )
  }
}


#' Compute summary statistics for a single group
#' @noRd
.compute_group_summary <- function(spectra, group_name, stats, na.rm = TRUE) {
  
  # Extract spectra for this group
  group_data <- get_spectra_by_group(spectra, group_name)
  spectra_list <- group_data$spectra
  
  # Get first spectrum to determine structure
  first_spectrum <- spectra_list[[1]]
  
  # Determine the spectral value column name
  spec_col <- .get_spectral_column(first_spectrum)
  
  # Get wavelengths
  wavelengths <- first_spectrum$wavelength
  n_wavelengths <- length(wavelengths)
  n_spectra <- length(spectra_list)
  
  # Start with wavelength column
  result <- data.frame(wavelength = wavelengths)
  
  # Handle single spectrum case
  if (n_spectra == 1) {
    values <- first_spectrum[[spec_col]]
    
    if ("mean" %in% stats) result$mean <- values
    if ("sd" %in% stats) result$sd <- rep(NA_real_, n_wavelengths)
    if ("se" %in% stats) result$se <- rep(NA_real_, n_wavelengths)
    if ("min" %in% stats) result$min <- values
    if ("max" %in% stats) result$max <- values
    if ("median" %in% stats) result$median <- values
    if ("n" %in% stats) result$n <- rep(1L, n_wavelengths)
    
    return(result)
  }
  
  # Extract spectral values into matrix (wavelengths × spectra)
  # Use do.call + cbind to ensure matrix output
  values_list <- lapply(spectra_list, function(s) s[[spec_col]])
  values_matrix <- do.call(cbind, values_list)
  
  # Compute requested statistics
  if ("mean" %in% stats) {
    result$mean <- rowMeans(values_matrix, na.rm = na.rm)
  }
  
  if ("sd" %in% stats) {
    result$sd <- apply(values_matrix, 1, sd, na.rm = na.rm)
  }
  
  if ("se" %in% stats) {
    sds <- apply(values_matrix, 1, sd, na.rm = na.rm)
    result$se <- sds / sqrt(n_spectra)
  }
  
  if ("min" %in% stats) {
    result$min <- apply(values_matrix, 1, min, na.rm = na.rm)
  }
  
  if ("max" %in% stats) {
    result$max <- apply(values_matrix, 1, max, na.rm = na.rm)
  }
  
  if ("median" %in% stats) {
    result$median <- apply(values_matrix, 1, median, na.rm = na.rm)
  }
  
  if ("n" %in% stats) {
    result$n <- rep(n_spectra, n_wavelengths)
  }
  
  return(result)
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


#' Create the output summary object
#' @noRd
.create_summary_object <- function(summary_list, group_names, stats,
                                   original_spectra) {
  
  # The summary list is the main content
  result <- summary_list
  
  # Set class
  class(result) <- c("ruffs_spectra_summary", "list")
  
  # Add metadata
  attr(result, "group_names") <- group_names
  attr(result, "stats_computed") <- stats
  
  # Preserve original metadata
  attr(result, "instrument") <- attr(original_spectra, "instrument")
  attr(result, "data_type") <- attr(original_spectra, "data_type")
  attr(result, "wavelength_range") <- attr(original_spectra, "wavelength_range")
  
  return(result)
}


#' Report summary results to console
#' @noRd
.report_summary <- function(group_names, stats) {
  
  n_groups <- length(group_names)
  
  message(
    "Computed statistics for ", n_groups, " group(s): ",
    paste(stats, collapse = ", ")
  )
}


#' Print Method for Spectra Summary
#'
#' @param x A `ruffs_spectra_summary` object.
#' @param ... Additional arguments (unused).
#'
#' @export
print.ruffs_spectra_summary <- function(x, ...) {
  
  group_names <- attr(x, "group_names")
  stats_computed <- attr(x, "stats_computed")
  n_groups <- length(group_names)
  
  # Header
  cat("RUFFS Spectral Summary\n")
  cat(strrep("-", 40), "\n")
  
  # Metadata
  cat("Number of groups: ", n_groups, "\n")
  cat("Statistics:       ", paste(stats_computed, collapse = ", "), "\n")
  
  instrument <- attr(x, "instrument")
  data_type <- attr(x, "data_type")
  if (!is.null(instrument)) cat("Instrument:       ", instrument, "\n")
  if (!is.null(data_type)) cat("Data type:        ", data_type, "\n")
  
  cat("\n")
  
  # Group details
  cat("Groups:\n")
  for (grp in group_names) {
    grp_data <- x[[grp]]
    wl_range <- range(grp_data$wavelength)
    n_points <- nrow(grp_data)
    
    # Get sample size if available
    if ("n" %in% names(grp_data)) {
      n_spectra <- grp_data$n[1]
      cat("  - '", grp, "': n = ", n_spectra,
          ", ", n_points, " wavelengths (",
          wl_range[1], "-", wl_range[2], " nm)\n", sep = "")
    } else {
      cat("  - '", grp, "': ", n_points, " wavelengths (",
          wl_range[1], "-", wl_range[2], " nm)\n", sep = "")
    }
  }
  
  cat("\nAccess with: summary[[\"group_name\"]] or summary$group_name\n")
  
  invisible(x)
}



#' Plot Spectral Summary with Error Bands
#'
#' @description
#' Creates a plot of mean spectra with optional error bands (SD or SE) for
#' one or more groups from a spectral summary.
#'
#' @param summary A `ruffs_spectra_summary` object created by
#'   [summarise_spectra()].
#' @param groups Character vector of group names to plot. If `NULL` (default),
#'   all groups are plotted.
#' @param error_type Character string specifying the error band type:
#'   `"se"` (default), `"sd"`, or `"none"`.
#' @param colours Character vector of colours for each group. If `NULL`,
#'   uses default colour palette.
#' @param alpha Numeric value (0-1) for error band transparency. Default is 0.2.
#' @param legend_position Position for legend. Options include `"topright"`,
#'   `"topleft"`, `"bottomright"`, `"bottomleft"`, or `"none"` to suppress.
#' @param xlim Numeric vector of length 2 specifying x-axis (wavelength) limits.
#'   If `NULL` (default), uses full wavelength range.
#' @param ylim Numeric vector of length 2 specifying y-axis limits.
#'   If `NULL` (default), automatically calculated to fit data and error bands.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the summary object.
#'
#' @examples
#' \dontrun{
#' summary <- summarise_spectra(grouped_spectra)
#'
#' # Plot all groups with SE bands
#' plot_summary(summary)
#'
#' # Plot specific groups with SD bands
#' plot_summary(summary,
#'   groups = c("site1_oak", "site2_oak"),
#'   error_type = "sd")
#'
#' # Plot without error bands
#' plot_summary(summary, error_type = "none")
#'
#' # Specify axis limits
#' plot_summary(summary, xlim = c(400, 900), ylim = c(0, 0.5))
#'
#' # Focus on red edge region
#' plot_summary(summary, xlim = c(680, 750), ylim = c(0.1, 0.6))
#'
#' # Full SWIR range
#' plot_summary(summary, xlim = c(1500, 2400))
#' }
#'
#' @export
plot_summary <- function(summary,
                         groups = NULL,
                         error_type = c("se", "sd", "none"),
                         colours = NULL,
                         alpha = 0.2,
                         legend_position = "topright",
                         xlim = NULL,
                         ylim = NULL,
                         ...) {
  
  # Validate inputs
  if (!inherits(summary, "ruffs_spectra_summary")) {
    stop("Input must be a 'ruffs_spectra_summary' object.", call. = FALSE)
  }
  
  error_type <- match.arg(error_type)
  
  # Get groups to plot
  all_groups <- attr(summary, "group_names")
  if (is.null(groups)) {
    groups <- all_groups
  } else {
    invalid <- setdiff(groups, all_groups)
    if (length(invalid) > 0) {
      stop(
        "Group(s) not found: ", paste(invalid, collapse = ", "), "\n",
        "Available groups: ", paste(all_groups, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  n_groups <- length(groups)
  
  # Set up colours
  if (is.null(colours)) {
    colours <- .default_palette(n_groups)
  } else if (length(colours) < n_groups) {
    colours <- rep_len(colours, n_groups)
  }
  
  # Check required columns exist
  stats_computed <- attr(summary, "stats_computed")
  if (!"mean" %in% stats_computed) {
    stop("Summary must include 'mean' statistic for plotting.", call. = FALSE)
  }
  if (error_type != "none" && !error_type %in% stats_computed) {
    warning(
      "Error type '", error_type, "' not available in summary. ",
      "Plotting without error bands.",
      call. = FALSE
    )
    error_type <- "none"
  }
  
  # Calculate y-axis limits if not provided
  if (is.null(ylim)) {
    ylim <- .calculate_plot_limits(summary, groups, error_type, xlim)
  }
  
  # Get x-axis (wavelength) from first group
  wavelengths <- summary[[groups[1]]]$wavelength
  
  # Set xlim to full range if not provided
  if (is.null(xlim)) {
    xlim <- range(wavelengths)
  }
  
  # Set up plot
  plot(wavelengths, summary[[groups[1]]]$mean,
       type = "n",
       xlab = "Wavelength (nm)",
       ylab = .get_ylab(summary),
       xlim = xlim,
       ylim = ylim,
       ...)
  
  # Plot each group
  for (i in seq_along(groups)) {
    grp <- groups[i]
    grp_data <- summary[[grp]]
    col <- colours[i]
    
    # Add error band if requested
    if (error_type != "none") {
      error_vals <- grp_data[[error_type]]
      .add_error_band(
        x = grp_data$wavelength,
        y = grp_data$mean,
        error = error_vals,
        col = col,
        alpha = alpha
      )
    }
    
    # Add mean line
    lines(grp_data$wavelength, grp_data$mean, col = col, lwd = 2)
  }
  
  # Add legend
  if (!identical(legend_position, "none")) {
    legend(legend_position,
           legend = groups,
           col = colours,
           lwd = 2,
           bty = "n")
  }
  
  invisible(summary)
}


#' Default colour palette
#' @noRd
.default_palette <- function(n) {
  if (n <= 8) {
    # Colourblind-friendly palette
    cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
              "#0072B2", "#D55E00", "#CC79A7", "#999999")
    return(cols[seq_len(n)])
  } else {
    return(rainbow(n))
  }
}


#' Calculate y-axis limits for plot
#' @noRd
.calculate_plot_limits <- function(summary, groups, error_type, xlim = NULL) {
  
  all_vals <- c()
  
  for (grp in groups) {
    grp_data <- summary[[grp]]
    
    # If xlim specified, subset data to that range for limit calculation
    if (!is.null(xlim)) {
      in_range <- grp_data$wavelength >= xlim[1] & grp_data$wavelength <= xlim[2]
      means <- grp_data$mean[in_range]
      if (error_type != "none" && error_type %in% names(grp_data)) {
        errors <- grp_data[[error_type]][in_range]
      } else {
        errors <- NULL
      }
    } else {
      means <- grp_data$mean
      if (error_type != "none" && error_type %in% names(grp_data)) {
        errors <- grp_data[[error_type]]
      } else {
        errors <- NULL
      }
    }
    
    if (!is.null(errors)) {
      all_vals <- c(all_vals, means - errors, means + errors)
    } else {
      all_vals <- c(all_vals, means)
    }
  }
  
  range(all_vals, na.rm = TRUE)
}


#' Get appropriate y-axis label
#' @noRd
.get_ylab <- function(summary) {
  data_type <- attr(summary, "data_type")
  
  if (is.null(data_type)) {
    return("Spectral Value")
  }
  
  switch(tolower(data_type),
         "reflectance" = "Reflectance",
         "radiance" = expression(paste("Radiance (W ", m^-2, " ", sr^-1, " ", nm^-1, ")")),
         "raw_dn" = "Digital Number (DN)",
         "Spectral Value"
  )
}


#' Add semi-transparent error band to plot
#' @noRd
.add_error_band <- function(x, y, error, col, alpha) {
  
  # Convert colour to RGB and add alpha
  rgb_col <- col2rgb(col) / 255
  band_col <- rgb(rgb_col[1], rgb_col[2], rgb_col[3], alpha = alpha)
  
  # Create polygon coordinates
  polygon(
    x = c(x, rev(x)),
    y = c(y + error, rev(y - error)),
    col = band_col,
    border = NA
  )
}