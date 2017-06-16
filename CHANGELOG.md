# Change Log
All notable changes to this project will be documented in this file. Project started
2016-04-01 but this log started 2017-03-16, after version 0.2.2. More details can be found
in the git logs.


## unreleased
### Changed 
### Removed
### Fixed
### Added


## [0.3.2] - 2017-06-16
### Changed 
### Removed
### Fixed
- reduced order of high-pass butterworth filter (with cut-off of 0.1Hz; see change in
  v0.3.1) from 5 to 2 in remove_artefacts.m function
### Added
- name of .mat file as input argument to 'resample_savemat.m'
- option in filter_butterworth_withnans.m to have separate filter orders for high and
  low-pass components of a bandpass filter
- file to plot rEEG (plot_rEEG.m)


## [0.3.1] - 2017-04-07
### Changed 
- bandpass filter changed from 0.5-40 Hz to 0.1-40 Hz before removing large-amplitude
  artefacts

### Removed
- bandpass filtering (0.5-40 Hz) on EEG prior to removing artefacts (constant values or
  sudden jumps in amplitude)

### Fixed
- avoid failing if only 1 channel for connectivity features (return NaN instead)
### Added


## [0.3.0] - 2017-03-29
### Changed 
- absolute and relative spectral power ('spectral\_power' and 'spectral\_relative\_power')
  now using periodogram (rectangular window) instead of Welch's power-spectral-density
  (PSD) estimate.
- all spectral estimates (PSD, robust-PSD, and periodogram) now used FFT length equal to
  window length (previously was in 2^n form, and at least >=256), i.e. no zero-padding.

### Removed
- PSD option for absolute and relative spectral power removed (replaced with periodogram).

### Fixed
- coherence function broken
- split of frequency bands for spectral and connectivity features when length of PSD is
  odd

### Added
- documentation to file headers
- this CHANGELOG.md!
- option to select PSD estimate for spectral features: either Welch PSD, robust-PSD
  (modified Welch PSD), or periodogram
- coherence function can use either PSD or robust-PSD 


## [0.2.2] - 2017-03-09
### Added
- surrogate data approach with coherence function
- NaN option to 'FILTER\_REPLACE\_ARTEFACTS'
- missing 'get\_window.m' and 'make\_odd.m' function
- now using semantic versioning (see http://semver.org/)

### Fixed
- scaling factor for robust-PSD

