#' Read SVC .sig files
#' 
#' Reads spectral data from Spectra Vista Corporation (SVC) ASCII .sig files.
#' 
#' @param file_path Character vector of file path(s) to .sig file(s). 
#'   If NULL, reads all .sig files from current working directory.
#' @param data_type Type of data to extract. Options:
#'   \describe{
#'     \item{"reflectance"}{Reflectance data (default)}
#'     \item{"radiance_target"}{Target radiance}
#'     \item{"radiance_reference"}{Reference radiance}
#'     \item{"irradiance"}{Irradiance data}
#'   }
#' @return A ruffs_spectra object (single file) or ruffs_spectra_collection (multiple files)
#' @export
#' @examples
#' \dontrun{
#' # Read single file
#' spec <- svcreader("sample.sig")
#' 
#' # Read all SIG files in current directory
#' specs <- svcreader()
#' 
#' # Read target radiance
#' rad <- svcreader("sample.sig", data_type = "radiance_target")
#' }
svcreader <- function(file_path = NULL, data_type = "reflectance") {
  
  # Validate data_type
  valid_types <- c("reflectance", "radiance_target", "radiance_reference", "irradiance")
  data_type <- tolower(data_type)
  if (!data_type %in% valid_types) {
    stop(sprintf("Invalid data_type '%s'. Must be one of: %s",
                 data_type, paste(valid_types, collapse = ", ")))
  }
  
  # Handle NULL file_path - use current directory
  if (is.null(file_path)) {
    file_path <- list.files(pattern = "\\.sig$", ignore.case = TRUE, 
                            full.names = TRUE)
    
    if (length(file_path) == 0) {
      stop("No .sig files found in current working directory")
    }
    
    message(sprintf("Found %d .sig file(s) in current directory", length(file_path)))
  }
  
  # Handle multiple files
  if (length(file_path) > 1) {
    spectra_list <- lapply(file_path, function(fp) {
      svcreader(fp, data_type = data_type)
    })
    names(spectra_list) <- basename(file_path)
    return(new_ruffs_spectra_collection(spectra_list))
  }
  
  # Single file processing
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  # Read file
  lines <- readLines(file_path, warn = FALSE)
  
  # Parse header
  header <- .svc_parse_header(lines)
  
  # Find and parse data section
  data_start <- grep("^data=", lines, ignore.case = TRUE)[1]
  if (is.na(data_start)) {
    stop("Could not find 'data=' line in SVC file")
  }
  
  data_lines <- lines[(data_start + 1):length(lines)]
  spectral_data <- .svc_parse_data(data_lines, data_type)
  
  # Determine column name
  col_name <- switch(data_type,
                     "reflectance" = "reflectance",
                     "radiance_target" = "radiance_target",
                     "radiance_reference" = "radiance_reference",
                     "irradiance" = "irradiance"
  )
  
  # Build metadata list
  metadata <- list(
    instrument_type = "SVC",
    filename = basename(file_path),
    filepath = normalizePath(file_path),
    data_type = data_type,
    instrument_model = header$instrument,
    instrument_number = header$inst_number,
    file_version = header$file_version,
    acquisition_time = .svc_parse_time(header$time),
    longitude = header$longitude,
    latitude = header$latitude,
    gpstime = header$gpstime,
    integration_time_vnir = header$it,
    integration_time_swir1 = header$it2,
    integration_time_swir2 = header$it3,
    temp_vnir = header$temp,
    temp_swir1 = header$temp2,
    temp_swir2 = header$temp3,
    battery_voltage = header$battery,
    measurement_type = header$comm,
    optic = header$optic,
    overlap_vnir_swir1 = c(990, 1010),
    overlap_swir1_swir2 = c(1790, 1810),
    memory_slot = header$memory_slot,
    factors = header$factors,
    comment = header$comment
  )
  
  # Create and return ruffs_spectra object
  spectra <- new_ruffs_spectra(
    wavelength = spectral_data$wavelength,
    values = spectral_data$values,
    value_name = col_name,
    metadata = metadata
  )
  
  return(spectra)
}

# Internal helper functions for SVC reading
.svc_parse_header <- function(lines) {
  header <- list()
  
  .extract <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)[1]
    if (length(line) == 0 || is.na(line)) return(NULL)
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) < 2) return(NULL)
    return(trimws(parts[2]))
  }
  
  header$instrument <- .extract("^name")
  header$inst_number <- .extract("^instrument")
  header$file_version <- .extract("^file version")
  header$time <- .extract("^time")
  header$longitude <- as.numeric(.extract("^longitude"))
  header$latitude <- as.numeric(.extract("^latitude"))
  header$gpstime <- .extract("^gpstime")
  header$it <- as.numeric(.extract("^integration="))
  
  # Extract all integration times from integration line
  int_line <- .extract("^integration=")
  if (!is.null(int_line)) {
    int_vals <- as.numeric(strsplit(int_line, ",")[[1]])
    if (length(int_vals) >= 3) {
      header$it <- int_vals[1]
      header$it2 <- int_vals[2]
      header$it3 <- int_vals[3]
    }
  }
  
  # Extract temperatures from temp line
  temp_line <- .extract("^temp=")
  if (!is.null(temp_line)) {
    temp_vals <- as.numeric(strsplit(temp_line, ",")[[1]])
    if (length(temp_vals) >= 3) {
      header$temp <- temp_vals[1]
      header$temp2 <- temp_vals[2]
      header$temp3 <- temp_vals[3]
    }
  }
  
  # Extract battery voltage (may have multiple values)
  battery_line <- .extract("^battery=")
  if (!is.null(battery_line)) {
    battery_vals <- as.numeric(strsplit(battery_line, ",")[[1]])
    header$battery <- battery_vals[1]
  }
  
  header$comm <- .extract("^comm")
  header$optic <- .extract("^optic")
  header$memory_slot <- .extract("^memory slot")
  header$comment <- .extract("^comment")
  
  # Parse factors
  factors_line <- .extract("^factors=")
  if (!is.null(factors_line)) {
    # Remove bracket content and extract numbers
    factors_clean <- gsub("\\[.*\\]", "", factors_line)
    header$factors <- as.numeric(strsplit(trimws(factors_clean), ",")[[1]])
  }
  
  return(header)
}

.svc_parse_data <- function(data_lines, data_type) {
  # Remove empty lines
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]
  
  # Parse data - handle both space and tab delimiters
  data_matrix <- do.call(rbind, lapply(data_lines, function(line) {
    # Split by whitespace (spaces or tabs)
    values <- strsplit(trimws(line), "\\s+")[[1]]
    as.numeric(values)
  }))
  
  # SVC column structure:
  # Column 1: wavelength
  # Column 2: target radiance
  # Column 3: reference radiance  
  # Column 4: reflectance (if present)
  
  wavelength <- data_matrix[, 1]
  
  values <- switch(data_type,
                   "reflectance" = {
                     if (ncol(data_matrix) >= 4) {
                       data_matrix[, 4] / 100
                     } else if (ncol(data_matrix) >= 3) {
                       # Calculate from target/reference
                       data_matrix[, 2] / data_matrix[, 3]
                     } else {
                       stop("Reflectance data not available in file")
                     }
                   },
                   "radiance_target" = {
                     if (ncol(data_matrix) >= 2) {
                       data_matrix[, 2]
                     } else {
                       stop("Target radiance not available in file")
                     }
                   },
                   "radiance_reference" = {
                     if (ncol(data_matrix) >= 3) {
                       data_matrix[, 3]
                     } else {
                       stop("Reference radiance not available in file")
                     }
                   },
                   "irradiance" = {
                     if (ncol(data_matrix) >= 5) {
                       data_matrix[, 5]
                     } else {
                       stop("Irradiance data not available in file")
                     }
                   }
  )
  
  return(list(wavelength = wavelength, values = values))
}

.svc_parse_time <- function(time_str) {
  if (is.null(time_str) || nchar(time_str) == 0) {
    return(NULL)
  }
  
  # SVC time format can have multiple timestamps separated by comma
  # Take the first one
  time_str <- strsplit(time_str, ",")[[1]][1]
  time_str <- trimws(time_str)
  
  tryCatch({
    time_obj <- as.POSIXct(time_str, format = "%m/%d/%Y %I:%M:%S%p", tz = "UTC")
    return(format(time_obj, "%Y-%m-%d %H:%M:%S"))
  }, error = function(e) {
    return(time_str)
  })
}