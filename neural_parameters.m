%---------------------------------------------------------------------
%% PREPROCESSING (lowpass filter and resample)
%---------------------------------------------------------------------
LP_fc = 30;  % low-pass filter cut-off
Fs_new = 64; % down-sample to Fs_new (Hz)


%---------------------------------------------------------------------
%% DIRECTORIES
%---------------------------------------------------------------------
% FILL IN AS APPROPRIATE:
EEG_DATA_DIR = '';
EEG_DATA_DIR_MATFILES = '';


%---------------------------------------------------------------------
%% MONTAGE
%---------------------------------------------------------------------
% bipolar montage for NICU babies:
BI_MONT = {{'F4','C4'},{'F3','C3'},{'C4','T4'},{'C3','T3'},{'C4','Cz'},{'Cz','C3'}, ...
         {'C4','O2'},{'C3','O1'}};

% or 
% $$$ BI_MONT = {{'F4','C4'},{'F3','C3'}};


%---------------------------------------------------------------------
%% ARTEFACTS
%---------------------------------------------------------------------
REMOVE_ART = 1; % simple proceedure to remove artefacts; 0 to turn off

% some default values used for preterm infants (<32 weeks of gestation)
ART_HIGH_VOLT   = 1500;   % in mirco Vs
ART_TIME_COLLAR = 10;     % time collar (in seconds) around high-amplitude artefact

ART_DIFF_VOLT = 200;          % in mirco Vs
ART_DIFF_TIME_COLLAR = 0.5;   % time collar (in seconds) around fast jumps
ART_DIFF_MIN_TIME = 0.1;      % min time (in seconds) for flat (continuous) trace to be artefact

ART_ELEC_CHECK = 1;   % minimum length required for electrode check (in seconds)

ART_REF_LOW_CORR = 0.15; % if mean correlation coefficent across referential channels  
                       % is < this value then remove


% what to replace artefacts with before filtering?
% options: 1) zeros ('zeros') 
%          2) linear interpolation ('linear_interp')
%          3) cubic spline interpolation ('cubic_interp')
%          4) NaN ('nans'): replace with cubic spline before filtering and then NaN's
%          after filtering
FILTER_REPLACE_ARTEFACTS = 'nans';

feat_params_st.amplitude.FILTER_REPLACE_ARTEFACTS = FILTER_REPLACE_ARTEFACTS;
feat_params_st.rEEG.FILTER_REPLACE_ARTEFACTS = FILTER_REPLACE_ARTEFACTS;
feat_params_st.connectivity.FILTER_REPLACE_ARTEFACTS = FILTER_REPLACE_ARTEFACTS;
feat_params_st.FD.FILTER_REPLACE_ARTEFACTS = FILTER_REPLACE_ARTEFACTS;

%---------------------------------------------------------------------
%% FEATURES
%---------------------------------------------------------------------
% list of all features:
all_features_list;

% or if just want a subset of features, set here:
% $$$ FEATURE_SET_ALL = { ...
% $$$     'spectral_relative_power' ...
% $$$     ,'spectral_entropy' ...
% $$$     ,'spectral_edge_frequency' ...
% $$$     ,'amplitude_total_power' ...
% $$$     ,'connectivity_coh_mean' ...
% $$$                 };


% band-pass filter in this band:
FREQ_BANDS = [0.5 4; 4 7; 7 13; 13 30]; 

% these bands often used for preterm infants (<32 weeks GA):
% $$$ FREQ_BANDS = [0.5 3; 3 8; 8 15; 15 30]; 

%---------------------------------------------------------------------
% A. spectral features
%---------------------------------------------------------------------
% how to estimate the spectrum for 'spectral_flatness', 'spectral_entropy',
% spectral_edge_frequency features:
% 1) PSD: estimate power spectral density (e.g. Welch periodogram)
% 2) robust-PSD: median (instead of mean) of spectrogram 
% 3) periodogram: magnitude of the discrete Fourier transform
feat_params_st.spectral.method = 'PSD'; 

% length of time-domain analysis window and overlap:
% (applies to 'spectral_power','spectral_relative_power',
%  'spectral_flatness', and 'spectral_diff' features)
feat_params_st.spectral.L_window = 2; % in seconds
feat_params_st.spectral.window_type = 'hamm'; % type of window
feat_params_st.spectral.overlap = 50; % overlap in percentage
feat_params_st.spectral.freq_bands = FREQ_BANDS;
feat_params_st.spectral.total_freq_bands = [FREQ_BANDS(1) FREQ_BANDS(end)];
feat_params_st.spectral.SEF = 0.95;  % spectral edge frequency



% fractal dimension (FD):
feat_params_st.FD.method = 'higuchi'; % method to estimate FD, either 'higuchi' or 'katz'
feat_params_st.FD.freq_bands = [FREQ_BANDS(1) FREQ_BANDS(end)];
% $$$ feat_params_st.FD.freq_bands = FREQ_BANDS;
feat_params_st.FD.qmax = 6;  % Higuchi method: max. value of k

%---------------------------------------------------------------------
% B. amplitude features
%---------------------------------------------------------------------
% $$$ feat_params_st.amplitude.freq_bands = [FREQ_BANDS(1) FREQ_BANDS(end)];
feat_params_st.amplitude.freq_bands = FREQ_BANDS;

% for rEEG (range-EEG, similar to aEEG) from [1]
%
% [1] D O’Reilly, MA Navakatikyan, M Filip, D Greene, & LJ Van Marter (2012). Peak-to-peak
% amplitude in neonatal brain monitoring of premature infants. Clinical Neurophysiology,
% 123(11), 2139–53.
%
% settings in [1]: window = 2 seconds; overlap = 0%; and no log-linear scale
feat_params_st.rEEG.L_window = 2; % in seconds
feat_params_st.rEEG.window_type = 'rect'; % type of window
feat_params_st.rEEG.overlap = 0; % overlap in percentage
feat_params_st.rEEG.APPLY_LOG_LINEAR_SCALE = 0; % use this scale (either 0 or 1)
feat_params_st.rEEG.freq_bands = FREQ_BANDS;


%---------------------------------------------------------------------
% C. connectivity features
%---------------------------------------------------------------------
% how to estimate the cross spectrum for the coherence function:
% 1) PSD: estimate power spectral density (e.g. Welch periodogram)
% 2) bartlett-PSD: Welch periodogram with 0% overlap and rectangular window
%    (necessary if using the analytic assessment of zero coherence, see below)
feat_params_st.connectivity.method = 'bartlett-PSD'; 

feat_params_st.connectivity.freq_bands = FREQ_BANDS;
feat_params_st.connectivity.L_window = 8; % PSD window in seconds
feat_params_st.connectivity.overlap = 75; % PSD window percentage overlap
feat_params_st.connectivity.window_type = 'hamm'; % PSD window type



% find lower coherence limit using either either a surrogate-data
% approach [1] or an analytic threshold [2]
% [1] Faes L, Pinna GD, Porta A, Maestri R, Nollo G (2004). Surrogate data analysis for
%     assessing the significance of the coherence function. IEEE Transactions on
%     Biomedical Engineering, 51(7):1156–1166.
% [2] Halliday, DM, Rosenberg, JR, Amjad, AM, Breeze, P, Conway, BA, &
%     Farmer, SF. (1995). A framework for the analysis of mixed time series/point
%     process data--theory and application to the study of physiological tremor, single
%     motor unit discharges and electromyograms. Progress in Biophysics and Molecular
%     Biology, 64(2–3), 237–278.
% 
% options for 'feat_params_st.connectivity.coherence_zero_level' are: 
% 1) 'surr' for [1]
% 2) 'analytic' for [2]
% 3) '' not to implement (no threshold)
feat_params_st.connectivity.coherence_zero_level = 'analytic';
% alpha value for null-hypothesis disribution cut-off:
feat_params_st.connectivity.coherence_zero_alpha = 0.05;
% number of iterations required to generate null-hypothesis distribution if 
% using surrogate data approach ([2]):
feat_params_st.connectivity.coherence_surr_iter = 500; 

%---------------------------------------------------------------------
%% SHORT-TIME ANALYSIS on EEG
%---------------------------------------------------------------------
EPOCH_LENGTH = 64;  % seconds
EPOCH_OVERLAP = 50; % percent

EPOCH_IGNORE_PRC_NANS = 50; % if epoch has ≥ EPOCH_IGNORE_PRC_NANS (percent) then ignore 
