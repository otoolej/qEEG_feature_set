# Change Log
All notable changes to this project will be documented in this file. Project started
2016-04-01 but this log started 2017-03-16, after version 0.2.2. More details can be found
in the git logs.

## unreleased
### Changed 
### Removed
### Fixed
### Added


## [0.4.4] - 26-08-2020
### Changed 
- Returns Table of feature values (per epoch and per channel) from function
  `generate_all_features` if input argument `return_feat_epoch=true`
### Removed
### Fixed
- Bug fix for `connectivity_features`; was: if channel labels (cell of strings) are empty
  then would assume (incorrectly) only 2 channels and would use index [1,2] for left,right
  hemisphere; now: if no channel labels throw warning
### Added
- Check in `spectral_features`: is `L_window` longer than signal (and if so throw
  warning)?


## [0.4.3] - 2020-04-22
### Changed 
- changed behaviour for `IBI_burst_prc` and `IBI_burst_number` features: now removing any
  NaNs before calculating this metric
### Removed
### Fixed
- bug in rand_phase.m when taking the inverse-FFT (random signal will time-reversed but
  with same initial x[0] value); won't affect the coherence values.
- for rectangular window in `filter_zerophase.m` don't divide by window length
- for IBI features fixed removal of minimal interval of 16 samples and set to 1/4 seconds
  instead
- avoid error if maximum/median IBI duration is zero
- bug fix for `spectral_diff` feature: extra zero was included in spectral difference measure
  before taking the median
### Added

## [0.4.2] - 2019-04-29
### Changed 
- connectivity analysis: remove channel pairs that have common electrode, e.g. 'C3-Cz' and
  'Cz-C4'
### Removed
### Fixed
- tested with Matlab R2019a
- fixed syntax problem for fir1 function for R2019a
### Added


## [0.4.1] - 2019-03-19
### Changed 
- minor tweak to the overlap-and-add analysis when generating features on an epoch of EEG;
  now include as much of the signal as possible and pad with NaNs when necessary.
### Removed
### Fixed
- tested with Matlab R2018b;
- bug when downsampling: always using the resample function even when not necessary;
- bug in PSD power feature if length-epoch < (2 * PSD-window-length); 
- bug when length-epoch = signal-length, was summing over frequency bands.
### Added


## [0.4.0] - 2019-01-09
### Changed 
- Syntax for specifying the zero-level threshold of the coherence function in
  neural_parameters.m. NB: this will break compatibility with previous versions.
### Removed
- 'robust-PSD' for coherence function (as does not satisfy the Schwartz inequality).
### Fixed
- Coherence function for surrogate analysis, which assumed that auto-PSDs Pxx and Pyy
  would not change when randomising the phase of the Fourier transform for the surrogate
  signals. This holds only when PSD is a periodogram. Assumption removed and now Pxx and
  Pyy generated for each iteration.
### Added
- Analytic zero-coherence threshold (Halliday et al. 1995) for coherence
  estimates. Estimates coherence using the Bartlett PSD.


## [0.3.4] - 2018-09-18
### Changed 
- Higuchi method for fractal dimension can now incorporate NaNs
### Removed
### Fixed
- tested with Matlab R2018a 
### Added

## [0.3.3] - 2017-11-16
### Changed 
### Removed
### Fixed
- tested with Matlab R2017a
### Added
- licence file
- DOI links from zenodo


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

