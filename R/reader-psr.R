#' Read Spectral Evolution .sed files
#' 
#' Reads spectral data from Spectral Evolution PSR series ASCII .sed files.
#' 
#' @param file_path Character vector of file path(s) to .sed file(s). 
#'   If NULL, reads all .sed files from current working directory.
#' @param data_type Type of data to extract. Options:
#'   \describe{
#'     \item{"reflectance"}{Reflectance data (default, converted from percentage)}
#'     \item{"radiance"}{Radiance data}
#'     \item{"irradiance"}{Irradiance data}
#'     \item{"transmittance"}{Transmittance data (converted from percentage)}
#'   }
#' @return A ruffs_spectra object (single file) or ruffs_spectra_collection (multiple files)
#' @export
#' @examples
#' \dontrun{
#' # Read single file
#' spec <- psrreader("sample.sed")
#' 
#' # Read all SED files in current directory
#' specs <- psrreader()
#' 
#' # Read transmittance data
#' trans <- psrreader("sample.sed", data_type = "transmittance")
#' }
psrreader <- function(file_path = NULL, data_type = "reflectance") {
  
  # Validate data_type
  valid_types <- c("reflectance", "radiance", "irradiance", "transmittance")
  data_type <- tolower(data_type)
  if (!data_type %in% valid_types) {
    stop(sprintf("Invalid data_type '%s'. Must be one of: %s",
                 data_type, paste(valid_types, collapse = ", ")))
  }
  
  # Handle NULL file_path - use current directory
  if (is.null(file_path)) {
    file_path <- list.files(pattern = "\\.sed$", ignore.case = TRUE, 
                            full.names = TRUE)
    
    if (length(file_path) == 0) {
      stop("No .sed files found in current working directory")
    }
    
    message(sprintf("Found %d .sed file(s) in current directory", length(file_path)))
  }
  
  # Handle multiple files
  if (length(file_path) > 1) {
    spectra_list <- lapply(file_path, function(fp) {
      psrreader(fp, data_type = data_type)
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
  header <- .psr_parse_header(lines)
  
  # Check data type compatibility
  if (!is.null(header$data_type) && tolower(header$data_type) != data_type) {
    warning(sprintf("File contains %s data but %s was requested",
                    header$data_type, data_type))
  }
  
  # Parse spectral data
  spectral_data <- .psr_parse_data(lines, data_type)
  
  # Determine column name
  col_name <- switch(data_type,
                     "reflectance" = "reflectance",
                     "radiance" = "radiance",
                     "irradiance" = "irradiance",
                     "transmittance" = "transmittance"
  )
  
  # Build metadata list
  metadata <- list(
    instrument_type = "PSR",
    filename = basename(file_path),
    filepath = normalizePath(file_path),
    data_type = data_type,
    instrument_model = header$instrument,
    serial_number = header$serial_number,
    acquisition_time = .psr_parse_time(header$timestamp),
    integration_time = header$integration_time,
    scans_to_average = header$scans,
    calibration_file = header$calibration,
    gps_latitude = header$latitude,
    gps_longitude = header$longitude,
    gps_altitude = header$altitude,
    foreoptic = header$foreoptic,
    dark_mode = header$dark_mode,
    reference_file = header$reference,
    temperature = header$temperature,
    data_units = header$units,
    original_data_type = header$data_type,
    comment = header$comment,
    operator = header$operator
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

# Internal helper functions for PSR reading
.psr_parse_header <- function(lines) {
  header <- list()
  
  .extract <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)[1]
    if (length(line) == 0 || is.na(line)) return(NULL)
    parts <- strsplit(line, ":")[[1]]
    if (length(parts) < 2) return(NULL)
    return(trimws(paste(parts[-1], collapse = ":")))
  }
  
  header$instrument <- .extract("^Instrument")
  header$serial_number <- .extract("^Serial")
  header$timestamp <- .extract("^Time")
  header$calibration <- .extract("^Calibrated Reference")
  header$foreoptic <- .extract("^Foreoptic")
  header$dark_mode <- .extract("^Dark Mode")
  header$reference <- .extract("^Reference")
  header$data_type <- .extract("^Measurement")
  header$units <- .extract("^Units")
  header$comment <- .extract("^Comment")
  header$operator <- .extract("^Operator")
  
  # Parse numeric fields
  int_time <- .extract("^Integration")
  if (!is.null(int_time)) {
    header$integration_time <- as.numeric(gsub("[^0-9.]", "", int_time))
  }
  
  scans <- .extract("^Averages")
  if (!is.null(scans)) {
    # Extract first number from comma-separated values
    scan_vals <- strsplit(scans, ",")[[1]]
    header$scans <- as.numeric(trimws(scan_vals[1]))
  }
  
  lat <- .extract("^Latitude")
  if (!is.null(lat) && lat != "n/a") {
    header$latitude <- as.numeric(gsub("[^0-9.-]", "", lat))
  } else {
    header$latitude <- NULL
  }
  
  lon <- .extract("^Longitude")
  if (!is.null(lon) && lon != "n/a") {
    header$longitude <- as.numeric(gsub("[^0-9.-]", "", lon))
  } else {
    header$longitude <- NULL
  }
  
  alt <- .extract("^Altitude")
  if (!is.null(alt) && alt != "n/a") {
    header$altitude <- as.numeric(gsub("[^0-9.-]", "", alt))
  } else {
    header$altitude <- NULL
  }
  
  temp <- .extract("^Temperature")
  if (!is.null(temp)) {
    # Extract first number from comma-separated values
    temp_vals <- strsplit(temp, ",")[[1]]
    header$temperature <- as.numeric(trimws(temp_vals[1]))
  }
  
  return(header)
}

.psr_parse_data <- function(lines, data_type) {
  # Find data section (line after "Data:")
  data_start <- grep("^Data:", lines, ignore.case = TRUE)
  
  if (length(data_start) == 0) {
    stop("Could not find data section in PSR file")
  }
  
  # Skip the header line (Chan.# Wvl ...) and start from actual data
  data_lines <- lines[(data_start + 2):length(lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]
  
  # Parse tab-delimited data
  data_matrix <- do.call(rbind, lapply(data_lines, function(line) {
    values <- strsplit(trimws(line), "\\s+")[[1]]
    as.numeric(values)
  }))

  wavelength <- data_matrix[, 2]  # Column 2 for wavelength
  
  values <- switch(data_type,
                   "reflectance" = {
                     if (ncol(data_matrix) >= 5) {
                       # Column 5 for reflectance, convert from percentage to decimal
                       refl_percent <- data_matrix[, 5]
                       if (max(refl_percent, na.rm = TRUE) > 1.5) {
                         refl_percent / 100
                       } else {
                         refl_percent
                       }
                     } else {
                       stop("Reflectance data not available in file")
                     }
                   },
                   "radiance" = {
                     if (ncol(data_matrix) >= 4) {
                       data_matrix[, 4]  # Target radiance
                     } else {
                       stop("Radiance data not available in file")
                     }
                   },
                   "irradiance" = {
                     if (ncol(data_matrix) >= 3) {
                       data_matrix[, 3]  # Reference irradiance
                     } else {
                       stop("Irradiance data not available in file")
                     }
                   },
                   "transmittance" = {
                     if (ncol(data_matrix) >= 5) {
                       trans_percent <- data_matrix[, 5]
                       if (max(trans_percent, na.rm = TRUE) > 1.5) {
                         trans_percent / 100
                       } else {
                         trans_percent
                       }
                     } else {
                       stop("Transmittance data not available in file")
                     }
                   }
  )
  
  return(list(wavelength = wavelength, values = values))
}

.psr_parse_time <- function(time_str) {
  if (is.null(time_str) || nchar(time_str) == 0) {
    return(NULL)
  }
  
  # PSR may have multiple timestamps separated by comma
  # Take the first one
  time_str <- strsplit(time_str, ",")[[1]][1]
  time_str <- trimws(time_str)
  
  # Try multiple time formats
  formats <- c(
    "%m/%d/%Y %H:%M:%S",
    "%Y-%m-%d %H:%M:%S",
    "%m-%d-%Y %H:%M:%S",
    "%d/%m/%Y %H:%M:%S"
  )
  
  for (fmt in formats) {
    tryCatch({
      time_obj <- as.POSIXct(time_str, format = fmt, tz = "UTC")
      if (!is.na(time_obj)) {
        return(format(time_obj, "%Y-%m-%d %H:%M:%S"))
      }
    }, error = function(e) NULL)
  }
  
  return(time_str)
}