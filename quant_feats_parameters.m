%---------------------------------------------------------------------
% PREPROCESSING (lowpass filter and resample)
%---------------------------------------------------------------------
LP_fc=30;  % low-pass filter cut-off
Fs_new=64; % ensure that (Fs/Fs_new) is integer


%---------------------------------------------------------------------
% DIRECTORIES
%---------------------------------------------------------------------
% FILL IN AS APPROPRIATE:
EEG_DATA_DIR='';
EEG_DATA_DIR='~/eeg_data/Rhodri_twin_study/';
EEG_DATA_DIR_MATFILES='';
EEG_DATA_DIR_MATFILES='~/eeg_data/Rhodri_twin_study/downsampled_matfiles/';


%---------------------------------------------------------------------
% MONTAGE
%---------------------------------------------------------------------
% bipolar montage for NICU babies:
BI_MONT={{'F4','C4'},{'F3','C3'},{'C4','T4'},{'C3','T3'},{'C4','Cz'},{'Cz','C3'}, ...
         {'C4','O2'},{'C3','O1'}};

% or 
% $$$ BI_MONT={{'F4','C4'},{'F3','C3'}};


%---------------------------------------------------------------------
% ARTEFACTS
%---------------------------------------------------------------------
REMOVE_ART=1; % simple proceedure to remove artefacts; 0 to turn off

% some default values used for preterm infants (<32 weeks of gestation)
ART_HIGH_VOLT  =1500;   % in mirco Vs
ART_TIME_COLLAR=10;     % time collar (in seconds) around high-amplitude artefact

ART_DIFF_VOLT=200;          % in mirco Vs
ART_DIFF_TIME_COLLAR=0.5;   % time collar (in seconds) around fast jumps
ART_DIFF_MIN_TIME=0.1;      % min time (in seconds) for flat (continuous) trace to be artefact

ART_ELEC_CHECK=1;   % minimum length required for electrode check (in seconds)



%---------------------------------------------------------------------
% FEATURES
%---------------------------------------------------------------------
% list of all features:
all_features_list;

% or if just want a subset of features, set here:
% $$$ FEATURE_SET_ALL={ ...
% $$$     'relative_spectral_power' ...
% $$$     ,'spectral_entropy' ...
% $$$     ,'spectral_edge_frequency' ...
% $$$     ,'amplitude_total_power' ...
% $$$     ,'connectivity_BSI' ...
% $$$     ,'connectivity_corr' ...    
% $$$                 };


% band-pass filter in this band:
% $$$ FREQ_BANDS=[0.5 4; 4 7; 7 13; 13 30]; 

% these bands often used for preterm infants (<32 weeks GA):
FREQ_BANDS=[0.1 3; 3 8; 8 15; 15 30]; 


% three types ways to generate spectrum:
% (applies to 'spectral_power' and 'relative_spectral_power' features)
% 1) PSD: estimate power spectral density (e.g. Welch periodgram)
% 2) spectogram: frequency marginal of spectrogram
% 3) med-spectogram: median (instead of mean) of spectrogram 
feat_params_st.spec.method='PSD'; 

% length of time-domain analysis window and overlap:
% (applies to 'spectral_power','relative_spectral_power',
%  'spectral_flatness', and 'spectral_diff' features)
feat_params_st.spec.L_window=8; % in seconds
feat_params_st.spec.window_type='hamm'; % type of window
feat_params_st.spec.overlap=75; % overlap in percentage
feat_params_st.spec.freq_bands=FREQ_BANDS;
feat_params_st.spec.total_freq_bands=[FREQ_BANDS(1) FREQ_BANDS(end)];
feat_params_st.spec.SEF=0.9;  % spectral edge frequency

% for amplitude features:
% $$$ feat_params_st.amplitude.freq_bands=[FREQ_BANDS(1) FREQ_BANDS(end)];
feat_params_st.amplitude.freq_bands=FREQ_BANDS;

% for connectivity features:
% feat_params_st.connectivity.freq_bands=[FREQ_BANDS(1) FREQ_BANDS(end)];
feat_params_st.connectivity.freq_bands=FREQ_BANDS;
feat_params_st.connectivity.PSD_window=2; % seconds
feat_params_st.connectivity.PSD_overlap=50; % seconds


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
feat_params_st.rEEG.APPLY_LOG_LINEAR_SCALE=1; % use this scale (either 0 or 1)
feat_params_st.rEEG.freq_bands=[1 20];
feat_params_st.rEEG.freq_bands=FREQ_BANDS;


%---------------------------------------------------------------------
% SHORT-TIME ANALYSIS on EEG
%---------------------------------------------------------------------
EPOCH_LENGTH=64;  % seconds
EPOCH_OVERLAP=50; % percent
