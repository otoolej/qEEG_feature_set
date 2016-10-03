%-------------------------------------------------------------------------------
% gen_test_EEGdata: generate EEG-like data as coloured Gaussian noise 
%
% Syntax: data_st=gen_test_EEGdata(dur,Fs)
%
% Inputs: 
%     dur:              duration of EEG-like data in seconds (default 300 seconds)
%     Fs:               sampling frequency
%     include_bipolar:  include the bipolar montage aswell as the referential  
%
% Outputs: 
%     data_st:          structure including EEG data (referential montage), sampling
%                       frequency, and channel labels
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 01-09-2016
%
% last update: Time-stamp: <2016-09-12 15:29:54 (otoolej)>
%-------------------------------------------------------------------------------
function data_st=gen_test_EEGdata(dur,Fs,include_bipolar)
if(nargin<1 || isempty(dur)), dur=60*2; end
if(nargin<2 || isempty(Fs)), Fs=256; end
if(nargin<3 || isempty(include_bipolar)), include_bipolar=0; end




N=floor(dur*Fs);


%---------------------------------------------------------------------
% generate white Gaussian noise and filter with moving average filter
% (referential montage has common signal)
%---------------------------------------------------------------------
L_ma=Fs/4;
L_ma_noise=ceil(Fs/16);
N_channels=9;
eeg_common=5.*filter(ones(1,L_ma)./L_ma,1,randn(1,N)); 

eeg_data_ref=filter(ones(1,L_ma)./L_ma,1,randn(N_channels,N)')' + ...
    (filter(ones(1,L_ma_noise)./L_ma_noise,1,randn(N_channels,N)')')./100 + ...
    repmat(eeg_common,[N_channels 1]);


data_st.eeg_data_ref=eeg_data_ref.*100;
data_st.Fs=Fs;


%---------------------------------------------------------------------
% generate bipolar and monopolar labels:
%---------------------------------------------------------------------
data_st.ch_labels_bi={{'F4','C4'},{'F3','C3'},{'C4','T4'},{'C3','T3'}, ...
                   {'C4','Cz'},{'Cz','C3'}, ...
                   {'C4','O2'},{'C3','O1'}};
data_st.ch_labels_ref=unique([data_st.ch_labels_bi{:}]);



%---------------------------------------------------------------------
% include bipolar montage as well?
%---------------------------------------------------------------------
if(include_bipolar)
    % generate bi-polar montage:
    [data_st.eeg_data,data_st.ch_labels]=set_bi_montage(data_st.eeg_data_ref, ...
                                                      data_st.ch_labels_ref, ...
                                                      data_st.ch_labels_bi);
end




% plot:
Dbplot=0;
if(Dbplot)
    
    %---------------------------------------------------------------------
    % EEG viewer installed?
    % (https://github.com/otoolej/eeg_viewer)
    %---------------------------------------------------------------------
    eeg_plotgui_withannos('signals',data_st.eeg_data_ref, ...
                          'Fs',data_st.Fs, ...
                          'channel_labels',data_st.ch_labels_ref, ...
                          'epoch_length',30);
end



