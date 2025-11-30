# RUFFS <img src="man/figures/logo.png" align="right" height="139" alt="" />

**R sUite For Field Spectroscopy** — A comprehensive toolkit for reading and processing spectral data from field spectrometers.

## Overview

RUFFS provides a unified workflow for working with spectral data from field spectrometers from the following manufacturers:

- **ASD** (Analytical Spectral Devices, Inc,) — `.asd` binary files
- **SVC** (Spectra Vista Corporation) — `.sig` files  
- **PSR** (Spectral Evolution) — `.sed` files

## Installation

You can install RUFFS from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("NERC-FSF/RUFFS")
```

Or using `pak`:

```r
# install.packages("pak")
pak::pak("NERC-FSF/RUFFS")
```

## Quick Start

```r
library(RUFFS)

# Read spectral files (auto-detects instrument type)
spectra <- reader("path/to/files")

# Or use instrument-specific readers
spectra <- asdreader("sample.asd")
spectra <- svcreader("sample.sig")
spectra <- psrreader("sample.sed")

# Read all files in a directory
spectra <- reader(instrument_type="asd/svc/psr", data_type = "reflectance/radiance")  # Reads all compatible files in working directory
```

## Processing Pipeline

RUFFS provides a complete spectral processing workflow:
## Reading Data

```r
# Read multiple files into a collection
specs <- reader(c("file1.asd", "file2.asd", "file3.asd"))

# Extract different data types
reflectance <- asdreader("sample.asd", data_type = "reflectance")
radiance <- asdreader("sample.asd", data_type = "radiance")
raw_dn <- asdreader("sample.asd", data_type = "raw_dn")
```

### Interpolation

Interpolate wavelength scales to user specified increments, e.g. 1 nm, rather than in default spectral interval of instrument:

```r
# Interpolate to 1nm spacing
spec_interp <- interpolation(spec, target_wavelengths = seq(350, 2500, 1))

# Apply to entire collection
specs_interp <- apply_to_collection(specs, interpolation, 
                                     target_wavelengths = seq(350, 2500, 1))
```

### Detector Overlap Stitching

Remove artefacts at the boundaries between detectors:

```r
spec_stitched <- overlap_stitching(spec)
```

### Jump Correction

Smooth step artefacts at detector transitions:

```r
spec_corrected <- jump_correct(spec)
```

### Grouping and Analysis

Group spectra by filename patterns and compute statistics:

```r
# Group by filename components
grouped <- group_spectra(specs, separator = "_", position = 2)

# Compute group statistics
summary <- summarise_spectra(grouped)

# Visualise with error bands
plot_summary(summary, error_type = "se")
```

## Complete Workflow Example

```r
library(RUFFS)

# Full processing pipeline
results <- reader("path/to/data") |>
  apply_to_collection(interpolation, target_wavelengths = seq(350, 2500, 1)) |>
  apply_to_collection(overlap_stitching) |>
  apply_to_collection(jump_correct) |>
  group_spectra(separator = "_", position = 2) |>
  summarise_spectra()

# Visualise
plot_summary(results, xlim = c(400, 900))
```

## Features

| Feature | Description |
|---------|-------------|
| **Multi-instrument support** | Read ASD, SVC, and PSR file formats |
| **Unified data structure** | S3 classes with preserved metadata |
| **Interpolation** | Linear, spline, and constant methods |
| **Overlap stitching** | Remove detector boundary artefacts |
| **Jump correction** | Smooth step discontinuities |
| **Grouping** | Flexible filename-based grouping |
| **Statistics** | Mean, SD, SE, min, max, median by group |
| **Visualisation** | Plot summaries with error bands |


## License

MIT © Robbie Ramsay
