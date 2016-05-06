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
% last update: Time-stamp: <2016-04-29 18:41:13 (otoolej)>
%-------------------------------------------------------------------------------
function [eeg_data,Fs]=resample_savemat(fname,channel_names)
if(nargin<1 || isempty(fname)), fname=[]; end
if(nargin<2 || isempty(channel_names)), channel_names=[]; end


DBplot=0;
DBplot_test=0;
SAVE_DATA=1;


quant_feats_parameters;

eeg_data=[];


%---------------------------------------------------------------------
% 1. load from EDF
%---------------------------------------------------------------------

% if multiple files, then concatenate:
if(iscell(fname))
    N_files=length(fname);
    for n=1:N_files
        fname{n}=strip_file_extension(fname{n});
        [data{n},ch_labels{n},Fs{n}]=edfread_nicolet([EEG_DATA_DIR fname{n} '.edf'], ...
                                                     channel_names);
        
        % check if labels and Fs are the same for all files:
        if(n==1)
            ch_labels_keep=ch_labels{n};
            Fs_keep=Fs{n};
        end
        if(~isequal(ch_labels{n},ch_labels_keep) || ~isequal(Fs{n},Fs_keep))
            error('different channel names or sampling frequency in .edf files');
        end
    end
    ch_labels=ch_labels_keep; Fs=Fs_keep;
    tends=cumsum(cellfun(@(x) size(x,2), data));
    data=cell2mat(data);

    if(DBplot_test)
        figure(9); clf; hold all;
        hx(1)=subplot(211); hold all;
        plot((1:size(data,2))./Fs,data');
        for k=1:N_files
            line([1 1].*tends(k)/Fs,ylim,'color','k');
        end
    end
else
    fname=strip_file_extension(fname);
    [data,ch_labels,Fs]=edfread_nicolet([EEG_DATA_DIR fname '.edf'],channel_names);
    
    if(DBplot_test)
        figure(9); clf; hold all;
        hx(1)=subplot(211); hold all;
        plot((1:size(data,2))./Fs,data');
    end    
    
    % SPECIAL CASE FOR THE MOBERG DATA: reverse polarity of C4
% $$$     ic4=find(strcmp(ch_labels,'C4'));
% $$$     data(ic4,:)=-data(ic4,:);
end


% convert from mono-polar to bi-polar montage:
data_ref=data; ch_labels_ref=ch_labels;
[data,ch_labels]=set_bi_montage(data,ch_labels,BI_MONT);
N_channels=size(data,1);


%---------------------------------------------------------------------
% 2. remove artefacts?
%---------------------------------------------------------------------
if(REMOVE_ART)
    data=remove_artefacts(data,ch_labels,Fs,data_ref,ch_labels_ref);
end



%---------------------------------------------------------------------
% 2. PRE-PROCESS (lowpass filter and resample)
%---------------------------------------------------------------------
% a. filter:
lp=LP_fc;
parfor n=1:N_channels
    if(~all(isnan(data(n,:))))
        d_mean=nanmean(data(n,:));
        s1_filt=filter_zerophase(data(n,:)-d_mean,Fs,lp,[],4001);
        data(n,:)=s1_filt+d_mean;
    end
end

% b. decimate:
if(isa(Fs/Fs_new,'integer'))
    idec=1:Fs/Fs_new:size(data,2);
    
    eeg_data=NaN(N_channels,length(idec));
    for n=1:N_channels
        eeg_data(n,:)=data(n,idec);
    end
    
else
    eeg_data=NaN(N_channels,ceil(size(data,2)*Fs_new/Fs));
    for n=1:N_channels
        eeg_data(n,:)=resample(data(n,:),Fs_new,Fs);
    end
end
Fs=Fs_new;


%---------------------------------------------------------------------
% 3. SAVE
%---------------------------------------------------------------------
if(SAVE_DATA)
    if(iscell(fname))
        fname_stub=num2str(cell2mat(regexp(fname{1},'[\d+]','match')));
        mfname=[EEG_DATA_DIR_MATFILES filesep fname_stub '.mat'];
    else
        mfname=[EEG_DATA_DIR_MATFILES filesep fname '.mat'];
    end
    
    time_now=now;
    save(mfname, 'eeg_data','ch_labels','Fs','time_now');
    
    if(DBplot_test)
        hx(2)=subplot(212); hold all;
        plot((1:size(eeg_data,2))./Fs,eeg_data');
        linkaxes(hx,'x');
    end
end




%---------------------------------------------------------------------
% 4. PLOT
%---------------------------------------------------------------------
if(DBplot)
  eeg_plotgui_withannos('signals',eeg_data,'Fs',Fs,'channel_labels',ch_labels,...
                        'bipolar_montage',-1,'epoch_length',60*10);
end


function fname=strip_file_extension(fname)
if(length(fname)>4 && strcmp(fname(end-3:end),'.edf'))
    fname=fname(1:end-4);
end
