#' Read ASD binary files
#' 
#' Reads spectral data from Analytical Spectral Devices (ASD) binary files (.asd).
#' 
#' @param file_path Character vector of file path(s) to .asd file(s). 
#'   If NULL, reads all .asd files from current working directory.
#' @param data_type Type of data to extract. Options:
#'   \describe{
#'     \item{"reflectance"}{Reflectance data (default, requires reference in file)}
#'     \item{"radiance"}{Radiance data}
#'     \item{"raw_dn"}{Raw digital numbers}
#'     \item{"reference"}{Reference spectrum}
#'   }
#' @return A ruffs_spectra object (single file) or ruffs_spectra_collection (multiple files)
#' @export
#' @examples
#' \dontrun{
#' # Read single file
#' spec <- asdreader("sample.asd")
#' 
#' # Read all ASD files in current directory
#' specs <- asdreader()
#' 
#' # Read radiance data
#' rad <- asdreader("sample.asd", data_type = "radiance")
#' 
#' # Read raw DN
#' raw <- asdreader("sample.asd", data_type = "raw_dn")
#' }
asdreader <- function(file_path = NULL, data_type = "reflectance") {
  
  # Validate data_type
  valid_types <- c("reflectance", "radiance", "raw_dn", "reference")
  data_type <- tolower(data_type)
  if (!data_type %in% valid_types) {
    stop(sprintf("Invalid data_type '%s'. Must be one of: %s",
                 data_type, paste(valid_types, collapse = ", ")))
  }
  
  # Handle NULL file_path - use current directory
  if (is.null(file_path)) {
    file_path <- list.files(pattern = "\\.asd$", ignore.case = TRUE, 
                            full.names = TRUE)
    
    if (length(file_path) == 0) {
      stop("No .asd files found in current working directory")
    }
    
    message(sprintf("Found %d .asd file(s) in current directory", length(file_path)))
  }
  
  # Handle multiple files
  if (length(file_path) > 1) {
    spectra_list <- lapply(file_path, function(fp) {
      asdreader(fp, data_type = data_type)
    })
    names(spectra_list) <- basename(file_path)
    return(new_ruffs_spectra_collection(spectra_list))
  }
  
  # Single file processing
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  # Open binary connection
  conn <- file(file_path, "rb")
  on.exit(close(conn))
  
  # Read header
  header <- .asd_read_header(conn)
  
  # Extract wavelengths
  wavelengths <- .asd_read_wavelengths(conn, header)
  
  # Read spectrum and white reference data
  spec_data <- .asd_read_spectrum_data(conn, header)
  
  # Process to requested data type
  data_values <- .asd_process_spectra(spec_data, header, data_type)
  
  # Determine column name
  col_name <- switch(data_type,
                     "reflectance" = "reflectance",
                     "radiance" = "radiance",
                     "raw_dn" = "raw_dn",
                     "reference" = "reference"
  )
  
  # Build metadata list
  metadata <- list(
    instrument_type = "ASD",
    filename = basename(file_path),
    filepath = normalizePath(file_path),
    data_type = data_type,
    file_version = header$file_version,
    file_size = file.info(file_path)$size,
    instrument_model = header$instrument,
    serial_number = header$serial_number,
    acquisition_time = .asd_parse_time(header$gps_time),
    integration_time = header$integration_time,
    foreoptic = header$foreoptic,
    dark_current_correction = header$dc_corr,
    measurement_type = header$data_type_name,
    data_format = header$data_format_name,
    channels = header$channels,
    reference_time = .asd_parse_time(header$reference_time),
    reference_description = header$white_reference,
    gps_latitude = header$gps_latitude,
    gps_longitude = header$gps_longitude,
    gps_altitude = header$gps_altitude,
    comment = header$comment,
    calibration_series = header$calibration_series,
    splice1_wavelength = header$splice1_wavelength,
    splice2_wavelength = header$splice2_wavelength,
    swir1_gain = header$swir1_gain,
    swir2_gain = header$swir2_gain
  )
  
  # Create and return ruffs_spectra object
  spectra <- new_ruffs_spectra(
    wavelength = wavelengths,
    values = data_values,
    value_name = col_name,
    metadata = metadata
  )
  
  return(spectra)
}

# Internal helper functions for ASD reading
.asd_read_header <- function(conn) {
  seek(conn, 0)
  header <- list()
  
  # Company name (3 bytes)
  seek(conn, 0)
  header$company <- readBin(conn, character(), n = 1, size = 3, endian = "little")
  
  # File version (1 byte at position 3)
  seek(conn, 3)
  header$file_version <- readBin(conn, "integer", n = 1, size = 1)
  
  # Comments (157 bytes at position 4)
  seek(conn, 4)
  header$comment <- readBin(conn, character(), n = 1, size = 157, endian = "little")
  
  # Program version (1 byte at position 178)
  seek(conn, 178)
  res <- readBin(conn, "integer", size = 1, n = 1)
  major <- bitwShiftR(res, 4)
  minor <- bitwAnd(res, 7)
  header$program_version <- paste0(major, '.', minor)
  
  # Format version (1 byte at position 179)
  seek(conn, 179)
  res <- readBin(conn, raw(), size = 1, n = 1)
  res <- as.integer(res)
  major <- bitwShiftR(res, 4)
  minor <- bitwAnd(res, 7)
  header$format_version <- paste0(major, '.', minor)
  
  # DC correction flag (1 byte at position 181)
  seek(conn, 181)
  header$dc_corr <- readBin(conn, "logical", n = 1, size = 1)
  
  # DC time (4 bytes at position 182)
  seek(conn, 182)
  header$dc_time <- readBin(conn, "integer", n = 1, size = 4, endian = "little")
  
  # Data type (1 byte at position 186)
  seek(conn, 186)
  data_type_code <- readBin(conn, "integer", n = 1, size = 1)
  data_types <- c('raw', 'reflectance', 'radiance', 'no_units', 'irradiance', 
                  'qi', 'transmittance', 'unknown', 'absorbance')
  header$data_type_name <- data_types[data_type_code + 1]
  
  # Reference time (4 bytes at position 187)
  seek(conn, 187)
  header$reference_time <- readBin(conn, "integer", n = 1, size = 4, endian = "little")
  
  # Starting wavelength (4 bytes at position 191)
  seek(conn, 191)
  header$ch1_wavel <- readBin(conn, "numeric", n = 1, size = 4, endian = "little")
  
  # Wavelength step (4 bytes at position 195)
  seek(conn, 195)
  header$wavel_step <- readBin(conn, "numeric", n = 1, size = 4, endian = "little")
  
  # Data format (1 byte at position 199)
  seek(conn, 199)
  data_format_code <- readBin(conn, "integer", n = 1, size = 1)
  data_formats <- c('numeric', 'integer', 'double', 'unknown')
  header$data_format <- data_format_code
  header$data_format_name <- data_formats[data_format_code + 1]
  
  # Channels (2 bytes at position 204)
  seek(conn, 204)
  header$channels <- readBin(conn, "integer", n = 1, size = 2, endian = "little")
  
  # GPS time (4 bytes at position 344)
  seek(conn, 344)
  header$gps_time <- readBin(conn, "integer", n = 1, size = 4, endian = "little")
  
  # Integration time (4 bytes at position 390)
  seek(conn, 390)
  header$integration_time <- readBin(conn, "integer", n = 1, size = 4, endian = "little")
  
  # Foreoptic (2 bytes at position 394)
  seek(conn, 394)
  header$foreoptic <- readBin(conn, "integer", n = 1, size = 2, endian = "little")
  
  # White reference description (32 bytes at position 356)
  seek(conn, 356)
  header$white_reference <- readBin(conn, character(), n = 1, size = 32, endian = "little")
  
  # Calibration series (2 bytes at position 398)
  seek(conn, 398)
  header$calibration_series <- readBin(conn, "integer", n = 1, size = 2, endian = "little")
  
  # Instrument number (2 bytes at position 400)
  seek(conn, 400)
  header$instrument_num <- readBin(conn, "integer", n = 1, size = 2, endian = "little")
  
  # Instrument type (1 byte at position 431)
  seek(conn, 431)
  inst_code <- readBin(conn, "integer", n = 1, size = 1)
  instruments <- c('unknown', 'PSII', 'LSVNIR', 'FieldSpec VNIR', 'FieldSpec FR', 
                   'FieldSpec NIR', 'CHEM', 'FieldSpec FullRange Unattended')
  header$instrument <- instruments[inst_code + 1]
  
  # Serial number (stored at position 396 as instrument_num for now)
  header$serial_number <- as.character(header$instrument_num)
  
  # GPS data flag (1 byte at position 433)
  seek(conn, 433)
  header$gps_data_flag <- readBin(conn, "integer", n = 1, size = 1)
  
  # GPS coordinates if available
  if (header$gps_data_flag == 1) {
    seek(conn, 334)
    gps_raw <- readBin(conn, raw(), n = 56, endian = "little")
    gps_con <- rawConnection(gps_raw)
    on.exit(close(gps_con), add = TRUE)
    
    header$true_heading <- readBin(gps_con, "numeric", 1, size = 8, endian = "little")
    header$speed <- readBin(gps_con, "numeric", 1, size = 8, endian = "little")
    header$gps_latitude <- readBin(gps_con, "numeric", 1, size = 8, endian = "little")
    header$gps_longitude <- readBin(gps_con, "numeric", 1, size = 8, endian = "little")
    header$gps_altitude <- readBin(gps_con, "numeric", 1, size = 8, endian = "little")
  } else {
    header$gps_latitude <- NULL
    header$gps_longitude <- NULL
    header$gps_altitude <- NULL
  }
  
  # SWIR gains (2 bytes each at 436, 438)
  seek(conn, 436)
  header$swir1_gain <- readBin(conn, "integer", n = 1, size = 2, endian = "little")
  seek(conn, 438)
  header$swir2_gain <- readBin(conn, "integer", n = 1, size = 2, endian = "little")
  
  # Splice wavelengths (4 bytes each at 444, 448)
  seek(conn, 444)
  header$splice1_wavelength <- readBin(conn, "numeric", n = 1, size = 4, endian = "little")
  seek(conn, 448)
  header$splice2_wavelength <- readBin(conn, "numeric", n = 1, size = 4, endian = "little")
  
  return(header)
}

.asd_read_wavelengths <- function(conn, header) {
  # Calculate wavelengths from start and step
  wavelengths <- header$ch1_wavel + (0:(header$channels - 1)) * header$wavel_step
  return(wavelengths)
}

.asd_read_spectrum_data <- function(conn, header) {
  # Determine data format size
  data_size <- switch(header$data_format_name,
                      "numeric" = 4,
                      "integer" = 4,
                      "double" = 8,
                      4)  # default
  
  # Seek to spectrum data (byte 484)
  seek(conn, where = 484)
  
  # Read spectrum based on format
  spectrum <- readBin(conn, what = header$data_format_name, 
                      n = header$channels, 
                      size = data_size,
                      endian = "little")
  
  # Calculate position for reference header
  ref_header_pos <- 484 + header$channels * data_size
  
  # Read reference header
  seek(conn, where = ref_header_pos)
  wr_flag <- readBin(conn, what = "integer", size = 2, endian = "little")
  
  # Read white reference time (8 bytes)
  wr_time <- readBin(conn, what = "integer", size = 8, endian = "little")
  
  # Read spectrum time (8 bytes)
  spec_time <- readBin(conn, what = "integer", size = 8, endian = "little")
  
  # Read spectrum description length (2 bytes)
  spec_description_length <- readBin(conn, what = "integer", size = 2, endian = "little")
  
  # Read spectrum description if present
  if (spec_description_length > 0) {
    spec_description <- readChar(conn, nchars = spec_description_length)
  } else {
    spec_description <- ""
  }
  
  # Read white reference
  wr <- readBin(conn, what = header$data_format_name, 
                n = header$channels,
                size = data_size,
                endian = "little")
  
  return(list(
    spectrum = spectrum,
    wr = wr,
    wr_flag = wr_flag,
    wr_time = wr_time,
    spec_time = spec_time,
    spec_description = spec_description
  ))
}

.asd_normalize_spectrum <- function(spec, header) {
  # Get wavelengths
  wl <- header$ch1_wavel + (0:(header$channels - 1)) * header$wavel_step
  
  # Find indices for each detector region
  idx1 <- which(wl <= header$splice1_wavelength)
  idx2 <- which(wl > header$splice1_wavelength & wl <= header$splice2_wavelength)
  idx3 <- which(wl > header$splice2_wavelength)
  
  # Normalize each region
  spec[idx1] <- spec[idx1] / header$integration_time
  spec[idx2] <- spec[idx2] * header$swir1_gain / 2048
  spec[idx3] <- spec[idx3] * header$swir2_gain / 2048
  
  return(spec)
}

.asd_process_spectra <- function(spec_data, header, data_type) {
  
  # If reflectance data requested
  if (data_type == 'reflectance') {
    if (header$data_type_name %in% c('radiance', 'reflectance')) {
      result <- spec_data$spectrum / spec_data$wr
    } else if (header$data_type_name == 'raw') {
      result <- .asd_normalize_spectrum(spec_data$spectrum, header) / 
        .asd_normalize_spectrum(spec_data$wr, header)
    } else {
      stop(paste0('File only contains data of type ', header$data_type_name, '.'))
    }
    
    # If radiance data requested
  } else if (data_type == 'radiance') {
    if (header$data_type_name %in% c('radiance', 'reflectance', 'irradiance')) {
      result <- spec_data$spectrum
    } else if (header$data_type_name == 'raw') {
      result <- .asd_normalize_spectrum(spec_data$spectrum, header)
    } else {
      stop(paste0('File only contains data of type ', header$data_type_name, '.'))
    }
    
    # If raw DN requested
  } else if (data_type == 'raw_dn') {
    if (header$data_type_name == 'raw') {
      result <- spec_data$spectrum
    } else {
      stop(paste0('File only contains data of type ', header$data_type_name, '. Cannot extract raw DN.'))
    }
    
    # If white reference requested
  } else if (data_type == 'reference') {
    if (header$data_type_name %in% c('radiance', 'reflectance', 'irradiance')) {
      result <- spec_data$wr
    } else if (header$data_type_name == 'raw') {
      result <- .asd_normalize_spectrum(spec_data$wr, header)
    } else {
      stop(paste0('Cannot extract white reference from ', header$data_type_name, '.'))
    }
    
  } else {
    stop('Invalid type of data requested.')
  }
  
  return(result)
}

.asd_parse_time <- function(time_int) {
  if (is.null(time_int) || time_int == 0) {
    return(NULL)
  }
  
  # ASD time is seconds since Unix epoch
  time_obj <- as.POSIXct(time_int, origin = "1970-01-01", tz = "UTC")
  return(format(time_obj, "%Y-%m-%d %H:%M:%S"))
}