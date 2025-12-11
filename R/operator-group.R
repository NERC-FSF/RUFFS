#' Group Spectra by Filename Pattern
#'
#' @description
#' Groups individual spectra within a RUFFS spectral collection based on
#' hierarchical naming patterns in filenames. This enables subsequent group-wise
#' operations such as averaging, calculating standard deviations, or other
#' statistical summaries within groups.
#'
#' @param spectra A `ruffs_spectra_collection` object containing multiple spectra.
#' @param separator Character string used as the delimiter in filenames.
#'   Default is `"_"`.
#' @param position Integer or integer vector indicating which occurrence(s) of 
#'   the separator to use for grouping. Spectra are grouped by the portion of 
#'   the filename to the LEFT of the specified separator occurrence(s). 
#'   Can be a single integer (e.g., `1`) or multiple integers (e.g., `c(1, 3)`).
#'   Default is `1`.
#' @param ignore_extension Logical. If `TRUE` (default), the file extension
#'   (e.g., `.sig`, `.sed`) is removed before parsing the filename.
#'
#' @return A `ruffs_spectra_collection_grouped` object, which is the input 
#'   collection with additional grouping information stored as attributes:
#'   \itemize{
#'     \item `group_names`: Character vector of unique group names
#'     \item `group_assignments`: Named character vector mapping each spectrum
#'       to its group
#'     \item `grouping_separator`: The separator used for grouping
#'     \item `grouping_position`: The position(s) used for grouping
#'   }
#'
#' @details
#' The function parses filenames to extract hierarchical group identifiers.
#' For example, with files named `1_tree_a.sig`, `1_tree_b.sig`, `2_grass_a.sig`:
#' \itemize{
#'   \item `separator = "_"`, `position = 1` groups by: `"1"`, `"1"`, `"2"`
#'   \item `separator = "_"`, `position = 2` groups by: `"1_tree"`, `"1_tree"`, `"2_grass"`
#'   \item `separator = "_"`, `position = c(1, 2)` groups by: `"1_tree"`, `"1_tree"`, `"2_grass"`
#' }
#' 
#' When multiple positions are specified, the function includes all text segments
#' between the first and last specified positions, including the separators.
#' For example, with `flower_a_R_`, using `position = c(1, 3)` groups by `"flower_R"`.
#'
#' Filenames that do not contain enough separator occurrences will trigger a
#' warning and be assigned to an `"ungrouped"` category.
#'
#' @examples
#' \dontrun{
#' # Read and process spectral collection
#' spectra <- reader("path/to/files")
#'
#' # Group by first underscore (e.g., "1" from "1_tree_a.sig")
#' grouped <- group_spectra(spectra, separator = "_", position = 1)
#'
#' # Group by second underscore (e.g., "1_tree" from "1_tree_a.sig")
#' grouped <- group_spectra(spectra, separator = "_", position = 2)
#'
#' # Group by first and third underscore (e.g., "flower_R" from "flower_a_R_")
#' grouped <- group_spectra(spectra, separator = "_", position = c(1, 3))
#'
#' # View group information
#' print(attr(grouped, "group_names"))
#' print(attr(grouped, "group_assignments"))
#' }
#'
#' @seealso [average_spectra()] for computing group means,
#'   [summarise_spectra()] for group-wise statistics
#'
#' @export
group_spectra <- function(spectra,
                          separator = "_",
                          position = 1,
                          ignore_extension = TRUE) {
  
  
  # Input validation
  .validate_group_inputs(spectra, separator, position)
  
  # Extract filenames from the collection
  
  filenames <- .get_spectra_names(spectra)
  
  if (length(filenames) == 0) {
    stop("No spectra found in the collection.", call. = FALSE)
  }
  
  # Remove file extensions if requested
  if (ignore_extension) {
    filenames_parsed <- .remove_extensions(filenames)
  } else {
    filenames_parsed <- filenames
  }
  
  # Extract group identifiers for each filename
  group_assignments <- .extract_group_names(
    filenames = filenames_parsed,
    original_names = filenames,
    separator = separator,
    position = position
  )
  
  # Get unique group names (preserving order of first occurrence)
  group_names <- unique(group_assignments)
  
  # Report grouping results
  .report_grouping(filenames, group_assignments, group_names)
  
  # Create grouped spectra object
  grouped_spectra <- .create_grouped_spectra(
    spectra = spectra,
    group_names = group_names,
    group_assignments = group_assignments,
    separator = separator,
    position = position
  )
  
  return(grouped_spectra)
}


#' Validate inputs for group_spectra
#' @noRd
.validate_group_inputs <- function(spectra, separator, position) {
  
  # Check spectra object - accept ruffs_spectra_collection
  if (!inherits(spectra, "ruffs_spectra_collection")) {
    stop(
      "Input must be a 'ruffs_spectra_collection' object.\n",
      "Received: ", class(spectra)[1],
      call. = FALSE
    )
  }
  
  # Check separator
  if (!is.character(separator) || length(separator) != 1) {
    stop(
      "'separator' must be a single character string.\n",
      "Received: ", typeof(separator), " of length ", length(separator),
      call. = FALSE
    )
  }
  
  if (nchar(separator) == 0) {
    stop("'separator' cannot be an empty string.", call. = FALSE)
  }
  
  # Check position - now allows vector of integers
  if (!is.numeric(position) || length(position) < 1) {
    stop(
      "'position' must be a single integer or vector of integers.\n",
      "Received: ", typeof(position), " of length ", length(position),
      call. = FALSE
    )
  }
  
  if (any(position < 1) || any(position != as.integer(position))) {
    stop(
      "'position' must be positive integer(s) (>= 1).\n",
      "Received: ", paste(position, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Check for duplicate positions
  if (length(position) > 1 && any(duplicated(position))) {
    stop(
      "'position' cannot contain duplicate values.\n",
      "Received: ", paste(position, collapse = ", "),
      call. = FALSE
    )
  }
}


#' Extract spectrum names from collection
#' @noRd
.get_spectra_names <- function(spectra) {
  
  # For ruffs_spectra_collection, the spectra are stored in a 'spectra' element
  # or as a named list directly
  if (!is.null(spectra$spectra)) {
    return(names(spectra$spectra))
  } else if (is.list(spectra)) {
    return(names(spectra))
  } else {
    # Fallback: try to get from attributes
    return(attr(spectra, "names"))
  }
}


#' Remove file extensions from filenames
#' @noRd
.remove_extensions <- function(filenames) {
  # Remove common spectral file extensions
  # Handles: .sig, .sed, .asd, .txt, etc.
  sub("\\.[^.]+$", "", filenames)
}


#' Extract group names based on separator position(s)
#' @noRd
.extract_group_names <- function(filenames, original_names, separator, position) {
  
  # Escape special regex characters in separator
  sep_escaped <- gsub("([.\\^$*+?{}\\[\\]|()])", "\\\\\\1", separator)
  
  # Sort positions to ensure we work from first to last
  position_sorted <- sort(position)
  
  group_assignments <- vapply(seq_along(filenames), function(i) {
    filename <- filenames[i]
    original <- original_names[i]
    
    # Find all positions of the separator
    matches <- gregexpr(sep_escaped, filename)[[1]]
    
    # Check if enough separators exist
    max_pos <- max(position_sorted)
    if (matches[1] == -1 || length(matches) < max_pos) {
      warning(
        "Filename '", original, "' has fewer than ", max_pos,
        " occurrence(s) of '", separator, "'. Assigned to 'ungrouped'.",
        call. = FALSE
      )
      return("ungrouped")
    }
    
    # Handle single or multiple positions
    if (length(position_sorted) == 1) {
      # Original behavior: extract to left of single position
      sep_position <- matches[position_sorted[1]]
      return(substr(filename, 1, sep_position - 1))
    } else {
      # Multiple positions: extract segments between first and last positions
      # and remove segments at excluded positions
      first_pos <- position_sorted[1]
      last_pos <- position_sorted[length(position_sorted)]
      
      # Get start (left of first position) and end (left of last position)
      start_idx <- matches[first_pos]
      end_idx <- matches[last_pos]
      
      # Extract the substring from start to end
      extracted <- substr(filename, 1, end_idx - 1)
      
      # Remove segments that should be excluded
      # Work backwards to preserve positions
      for (pos_idx in seq_along(matches)) {
        if (pos_idx > first_pos && pos_idx < last_pos && !(pos_idx %in% position_sorted)) {
          # This separator position should be removed along with text before it
          # Find the previous kept separator
          prev_kept <- max(position_sorted[position_sorted < pos_idx])
          next_kept <- min(position_sorted[position_sorted > pos_idx])
          
          # Remove from previous separator to this one
          prev_match <- matches[prev_kept]
          curr_match <- matches[pos_idx]
          
          # Remove the segment between prev_kept and pos_idx
          before <- substr(extracted, 1, prev_match - 1)
          after <- substr(extracted, curr_match, nchar(extracted))
          extracted <- paste0(before, after)
        }
      }
      
      return(extracted)
    }
    
  }, character(1))
  
  # Name the assignments with original filenames
  names(group_assignments) <- original_names
  
  return(group_assignments)
}


#' Report grouping results to console
#' @noRd
.report_grouping <- function(filenames, group_assignments, group_names) {
  n_spectra <- length(filenames)
  n_groups <- length(group_names)
  
  # Count spectra per group
  group_counts <- table(group_assignments)
  
  message(
    "Grouped ", n_spectra, " spectra into ", n_groups, " group(s):"
  )
  
  # Display group summary
  for (grp in group_names) {
    count <- group_counts[grp]
    message("  - '", grp, "': ", count, " spectrum/spectra")
  }
}


#' Create the grouped spectra object
#' @noRd
.create_grouped_spectra <- function(spectra, group_names, group_assignments,
                                    separator, position) {
  
  # Preserve existing attributes
  existing_attrs <- attributes(spectra)
  
  # Add grouping attributes
  attr(spectra, "group_names") <- group_names
  attr(spectra, "group_assignments") <- group_assignments
  attr(spectra, "grouping_separator") <- separator
  attr(spectra, "grouping_position") <- position
  attr(spectra, "grouped") <- TRUE
  
  # Add the grouped class while preserving existing classes
  class(spectra) <- c("ruffs_spectra_collection_grouped", class(spectra))
  
  return(spectra)
}

#' Get Group Names from Grouped Spectra
#'
#' @description
#' Retrieves the unique group names from a grouped spectral collection.
#'
#' @param spectra A `ruffs_spectra_collection_grouped` object.
#'
#' @return Character vector of group names.
#'
#' @export
get_group_names <- function(spectra) {
  if (!inherits(spectra, "ruffs_spectra_collection_grouped")) {
    stop("Input must be a grouped spectra collection. Use group_spectra() first.",
         call. = FALSE)
  }
  attr(spectra, "group_names")
}


#' Get Group Assignments from Grouped Spectra
#'
#' @description
#' Retrieves the mapping of spectrum names to group assignments.
#'
#' @param spectra A `ruffs_spectra_collection_grouped` object.
#'
#' @return Named character vector where names are spectrum filenames and
#'   values are group assignments.
#'
#' @export
get_group_assignments <- function(spectra) {
  if (!inherits(spectra, "ruffs_spectra_collection_grouped")) {
    stop("Input must be a grouped spectra collection. Use group_spectra() first.",
         call. = FALSE)
  }
  attr(spectra, "group_assignments")
}


#' Get Spectra by Group
#'
#' @description
#' Extracts all spectra belonging to a specific group.
#'
#' @param spectra A `ruffs_spectra_collection_grouped` object.
#' @param group_name Character string specifying which group to extract.
#'
#' @return A `ruffs_spectra_collection` object containing only spectra from the
#'   specified group.
#'
#' @export
get_spectra_by_group <- function(spectra, group_name) {
  if (!inherits(spectra, "ruffs_spectra_collection_grouped")) {
    stop("Input must be a grouped spectra collection. Use group_spectra() first.",
         call. = FALSE)
  }
  
  group_assignments <- attr(spectra, "group_assignments")
  group_names <- attr(spectra, "group_names")
  
  if (!group_name %in% group_names) {
    stop(
      "Group '", group_name, "' not found.\n",
      "Available groups: ", paste(group_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Get names of spectra in this group
  members <- names(group_assignments)[group_assignments == group_name]
  
  # Extract subset from the collection
  # Handle both $spectra and direct list structures
  if (!is.null(spectra$spectra)) {
    subset_spectra <- spectra$spectra[members]
  } else {
    subset_spectra <- spectra[members]
  }
  
  # Create new collection object
  result <- list(spectra = subset_spectra)
  class(result) <- "ruffs_spectra_collection"
  
  # Preserve key attributes
  attr(result, "instrument") <- attr(spectra, "instrument")
  attr(result, "data_type") <- attr(spectra, "data_type")
  attr(result, "wavelength_range") <- attr(spectra, "wavelength_range")
  
  return(result)
}


#' List Spectra in Each Group
#'
#' @description
#' Returns a list showing which spectra belong to each group.
#'
#' @param spectra A `ruffs_spectra_collection_grouped` object.
#'
#' @return A named list where names are group names and values are
#'   character vectors of spectrum filenames.
#'
#' @export
list_groups <- function(spectra) {
  if (!inherits(spectra, "ruffs_spectra_collection_grouped")) {
    stop("Input must be a grouped spectra collection. Use group_spectra() first.",
         call. = FALSE)
  }
  
  group_assignments <- attr(spectra, "group_assignments")
  group_names <- attr(spectra, "group_names")
  
  # Create list of members for each group
  group_list <- lapply(group_names, function(grp) {
    names(group_assignments)[group_assignments == grp]
  })
  names(group_list) <- group_names
  
  return(group_list)
}


#' Print Method for Grouped Spectra Collection
#'
#' @param x A `ruffs_spectra_collection_grouped` object.
#' @param ... Additional arguments (unused).
#'
#' @export
print.ruffs_spectra_collection_grouped <- function(x, ...) {
  
  # Get grouping info
  group_names <- attr(x, "group_names")
  group_assignments <- attr(x, "group_assignments")
  separator <- attr(x, "grouping_separator")
  position <- attr(x, "grouping_position")
  
  n_spectra <- length(group_assignments)
  n_groups <- length(group_names)
  
  # Header
  cat("RUFFS Grouped Spectral Collection\n")
  cat(strrep("-", 40), "\n")
  
  # Collection info
  cat("Total spectra:    ", n_spectra, "\n")
  cat("Number of groups: ", n_groups, "\n")
  
  # Format position display
  if (length(position) == 1) {
    pos_text <- paste0("#", position)
  } else {
    pos_text <- paste0("#", paste(position, collapse = ", #"))
  }
  cat("Grouping rule:    Left of '", separator, "' occurrence ", pos_text, "\n",
      sep = "")
  
  # Instrument and data type if available
  instrument <- attr(x, "instrument")
  data_type <- attr(x, "data_type")
  if (!is.null(instrument)) cat("Instrument:       ", instrument, "\n")
  if (!is.null(data_type)) cat("Data type:        ", data_type, "\n")
  
  cat("\n")
  
  # Group summary
  cat("Groups:\n")
  group_counts <- table(group_assignments)[group_names]  # Preserve order
  
  for (i in seq_along(group_names)) {
    grp <- group_names[i]
    count <- group_counts[i]
    members <- names(group_assignments)[group_assignments == grp]
    
    cat("  [", i, "] '", grp, "' (n=", count, ")\n", sep = "")
    
    # Show first few members
    n_show <- min(3, length(members))
    for (j in seq_len(n_show)) {
      cat("      - ", members[j], "\n", sep = "")
    }
    if (length(members) > n_show) {
      cat("      ... and ", length(members) - n_show, " more\n", sep = "")
    }
  }
  
  invisible(x)
}
