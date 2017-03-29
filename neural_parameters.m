%---------------------------------------------------------------------
%% PREPROCESSING (lowpass filter and resample)
%---------------------------------------------------------------------
LP_fc=30;  % low-pass filter cut-off
Fs_new=64; % down-sample to Fs_new (Hz)


%---------------------------------------------------------------------
%% DIRECTORIES
%---------------------------------------------------------------------
% FILL IN AS APPROPRIATE:
EEG_DATA_DIR='';
EEG_DATA_DIR_MATFILES='';


%---------------------------------------------------------------------
%% MONTAGE
%---------------------------------------------------------------------
% bipolar montage for NICU babies:
BI_MONT={{'F4','C4'},{'F3','C3'},{'C4','T4'},{'C3','T3'},{'C4','Cz'},{'Cz','C3'}, ...
         {'C4','O2'},{'C3','O1'}};

% or 
% $$$ BI_MONT={{'F4','C4'},{'F3','C3'}};


%---------------------------------------------------------------------
%% ARTEFACTS
%---------------------------------------------------------------------
REMOVE_ART=1; % simple proceedure to remove artefacts; 0 to turn off

% some default values used for preterm infants (<32 weeks of gestation)
ART_HIGH_VOLT  =1500;   % in mirco Vs
ART_TIME_COLLAR=10;     % time collar (in seconds) around high-amplitude artefact

ART_DIFF_VOLT=200;          % in mirco Vs
ART_DIFF_TIME_COLLAR=0.5;   % time collar (in seconds) around fast jumps
ART_DIFF_MIN_TIME=0.1;      % min time (in seconds) for flat (continuous) trace to be artefact

ART_ELEC_CHECK=1;   % minimum length required for electrode check (in seconds)

ART_REF_LOW_CORR=0.15; % if mean correlation coefficent across referential channels  
                       % is < this value then remove


% what to replace artefacts with before filtering?
% options: 1) zeros ('zeros') 
%          2) linear interpolation ('linear_interp')
%          3) cubic spline interpolation ('cubic_interp')
%          4) NaN ('nans'): replace with cubic spline before filtering and then NaN's
%          after filtering
FILTER_REPLACE_ARTEFACTS='cubic_interp';

feat_params_st.amplitude.FILTER_REPLACE_ARTEFACTS=FILTER_REPLACE_ARTEFACTS;
feat_params_st.rEEG.FILTER_REPLACE_ARTEFACTS=FILTER_REPLACE_ARTEFACTS;
feat_params_st.connectivity.FILTER_REPLACE_ARTEFACTS=FILTER_REPLACE_ARTEFACTS;
feat_params_st.FD.FILTER_REPLACE_ARTEFACTS=FILTER_REPLACE_ARTEFACTS;


%---------------------------------------------------------------------
%% FEATURES
%---------------------------------------------------------------------
% list of all features:
all_features_list;

% or if just want a subset of features, set here:
% $$$ FEATURE_SET_ALL={ ...
% $$$     'spectral_relative_power' ...
% $$$     ,'spectral_entropy' ...
% $$$     ,'spectral_edge_frequency' ...
% $$$     ,'amplitude_total_power' ...
% $$$     ,'connectivity_BSI' ...
% $$$     ,'connectivity_corr' ...    
% $$$                 };


% band-pass filter in this band:
FREQ_BANDS=[0.5 4; 4 7; 7 13; 13 30]; 

% these bands often used for preterm infants (<32 weeks GA):
% $$$ FREQ_BANDS=[0.5 3; 3 8; 8 15; 15 30]; 

%---------------------------------------------------------------------
% A. spectral features
%---------------------------------------------------------------------
% how to estimate the spectrum for 'spectral_flatness', 'spectral_entropy',
% spectral_edge_frequency features:
% 1) PSD: estimate power spectral density (e.g. Welch periodgram)
% 2) robust-PSD: median (instead of mean) of spectrogram 
% 3) periodogram: magnitude of the discrete Fourier transform
feat_params_st.spectral.method='PSD'; 

% length of time-domain analysis window and overlap:
% (applies to 'spectral_power','spectral_relative_power',
%  'spectral_flatness', and 'spectral_diff' features)
feat_params_st.spectral.L_window=2; % in seconds
feat_params_st.spectral.window_type='hamm'; % type of window
feat_params_st.spectral.overlap=50; % overlap in percentage
feat_params_st.spectral.freq_bands=FREQ_BANDS;
feat_params_st.spectral.total_freq_bands=[FREQ_BANDS(1) FREQ_BANDS(end)];
feat_params_st.spectral.SEF=0.95;  % spectral edge frequency



% fractal dimension (FD):
feat_params_st.FD.method='higuchi'; % method to estimate FD, either 'higuchi' or 'katz'
feat_params_st.FD.freq_bands=[FREQ_BANDS(1) FREQ_BANDS(end)];
% $$$ feat_params_st.FD.freq_bands=FREQ_BANDS;
feat_params_st.FD.qmax=6;  % Higuchi method: max. value of k

%---------------------------------------------------------------------
% B. amplitude features
%---------------------------------------------------------------------
% $$$ feat_params_st.amplitude.freq_bands=[FREQ_BANDS(1) FREQ_BANDS(end)];
feat_params_st.amplitude.freq_bands=FREQ_BANDS;

% for rEEG (range-EEG, similar to aEEG) from [1]
%
% [1] D O’Reilly, MA Navakatikyan, M Filip, D Greene, & LJ Van Marter (2012). Peak-to-peak
% amplitude in neonatal brain monitoring of premature infants. Clinical Neurophysiology,
% 123(11), 2139–53.
%
% settings in [1]: window = 2 seconds; overlap = 0%; and no log-linear scale
feat_params_st.rEEG.L_window=2; % in seconds
feat_params_st.rEEG.window_type='rect'; % type of window
feat_params_st.rEEG.overlap=0; % overlap in percentage
feat_params_st.rEEG.APPLY_LOG_LINEAR_SCALE=0; % use this scale (either 0 or 1)
feat_params_st.rEEG.freq_bands=FREQ_BANDS;


%---------------------------------------------------------------------
% C. connectivity features
%---------------------------------------------------------------------
% how to estimate the cross spectrum for the coherence function:
% 1) PSD: estimate power spectral density (e.g. Welch periodgram)
% 2) robust-PSD: median (instead of mean) of spectrogram 
feat_params_st.connectivity.method='PSD'; 

feat_params_st.connectivity.freq_bands=FREQ_BANDS;
feat_params_st.connectivity.L_window=2; % PSD window in seconds
feat_params_st.connectivity.overlap=50; % PSD window percentage overlap
feat_params_st.connectivity.window_type='hamm'; % PSD window type
% find lower coherence limit using surrogate data?
% (number of iterations required to generate null-hypothesis distribution;
%  set to 0 to turn off)
feat_params_st.connectivity.coherence_surr_data=100; 
% alpha value for null-hypothesis disribution cut-off:
feat_params_st.connectivity.coherence_surr_alpha=0.05;


%---------------------------------------------------------------------
%% SHORT-TIME ANALYSIS on EEG
%---------------------------------------------------------------------
EPOCH_LENGTH=64;  % seconds
EPOCH_OVERLAP=50; % percent

EPOCH_IGNORE_PRC_NANS=50; % if epoch has ≥ EPOCH_IGNORE_PRC_NANS (percent) then ignore 
