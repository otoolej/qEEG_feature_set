%-------------------------------------------------------------------------------
% view_eegdata: a) read in EEG from .edf files
%               b) remove artefacts
%               c) band-pass filter
%               d) downsample
%               e) save as .mat file
%
%
% Syntax: []=view_eegdata(fname,data,Fs)
%
% Inputs: 
%     fname,data,Fs - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%
% REQUIRES: edfread_nicolet.m; filter_zerophase.m

% John M. O' Toole, University College Cork
% Started: 27-05-2013
%
% last update: Time-stamp: <2016-04-07 13:39:05 (otoolej)>
%-------------------------------------------------------------------------------
function [data,Fs]=resample_savemat(fname,channel_names)
if(nargin<1 || isempty(fname)), fname=[]; end
if(nargin<2 || isempty(channel_names)), channel_names=[]; end


DBplot=0;
SAVE_DATA=1;


quant_feats_parameters;

eeg_data=[];


%---------------------------------------------------------------------
% 1. load from EDF
%---------------------------------------------------------------------
if(length(fname)>4 && ~isempty(find(ismember({'.mat','.edf'}, fname(end-3:end)))))
    fname=fname(1:end-4);
end

[data,ch_labels,Fs]=edfread_nicolet([EEG_DATA_DIR fname '.edf'],channel_names);

% convert from mono-polar to bi-polar montage:
[data,ch_labels]=set_bi_montage(data,ch_labels,BI_MONT);
N_channels=size(data,1);


%---------------------------------------------------------------------
% 2. remove artefacts?
%---------------------------------------------------------------------
if(REMOVE_ART)
    data=remove_artefacts(data,ch_labels,Fs);
end



%---------------------------------------------------------------------
% 2. PRE-PROCESS (lowpass filter and resample)
%---------------------------------------------------------------------
% a. filter:
d_mean=mean(data(1,:));
s1_filt=filter_zerophase(data(1,:),Fs,LP_fc,[],4001);
    
data_filt=zeros(N_channels,length(s1_filt));
data_filt(1,:)=s1_filt+d_mean;
for i=2:N_channels
    d_mean=mean(data(i,:));
    s1_filt=filter_zerophase(data(i,:),Fs,LP_fc,[],4001);
    data_filt(i,:)=s1_filt+d_mean;
end

% b. decimate:
eeg_data=data_filt(:,1:Fs/Fs_new:end);
Fs=Fs_new;


%---------------------------------------------------------------------
% 3. SAVE
%---------------------------------------------------------------------
if(SAVE_DATA)
    time_now=now;
    save([EEG_DATA_DIR_MATFILES filesep fname '.mat'], ...
         'eeg_data','ch_labels','Fs','time_now');
end




%---------------------------------------------------------------------
% 4. PLOT
%---------------------------------------------------------------------
if(DBplot)
  eeg_plotgui_withannos('signals',eeg_data,'Fs',Fs,'channel_labels',ch_labels,...
                        'bipolar_montage',-1,'epoch_length',60*10);
end
