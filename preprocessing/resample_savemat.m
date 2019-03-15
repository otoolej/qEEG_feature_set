%-------------------------------------------------------------------------------
% resample_savemat: a) read in EEG from .edf files
%                   b) remove artefacts
%                   c) band-pass filter
%                   d) downsample
%                   e) save as .mat file
%
%
% Syntax: [eeg_data,Fs]=resample_savemat(fname,channel_names)
%
% Inputs: 
%     fname         - EDF file name
%     channel_names - cell of channel names to read in; 
%
% Outputs: 
%     eeg_data - EEG data in bipolar montage
%     Fs       - sample frequency (in Hz)
%
%
% REQUIRES: user-supplied function to read in raw EEG

% John M. O' Toole, University College Cork
% Started: 27-05-2013
%
% last update: Time-stamp: <2019-03-15 10:17:52 (otoolej)>
%-------------------------------------------------------------------------------
function [eeg_data,Fs]=resample_savemat(fname,channel_names,fname_out)
if(nargin<1 || isempty(fname)), fname=[]; end
if(nargin<2 || isempty(channel_names)), channel_names=[]; end
if(nargin<3 || isempty(fname_out)), fname_out=[]; end



DBplot=0;
% EEG_viewer required for DBplot=1 
% (https://github.com/otoolej/eeg_viewer.git)
DBplot_test=0;
SAVE_DATA=1;


neural_parameters;

eeg_data=[]; Fs=[];


%---------------------------------------------------------------------
% 0. supply own function to read in raw EEG files (e.g. EDF files)
%
% input arguments:  file_name, channel_names
% output arguments: data (EEG data matrix), ch_labels (channel labels), 
%                   Fs (sampling frequency)
% 
%---------------------------------------------------------------------
% a) COMMENT:
fprintf('must supply own function to read in raw EEG; see %s (line 51)\n', ...
        mfilename); %('fullpath'));
return;
% b) UNCOMMENT and REPLACE <edfread_nicolet> with own function to read in raw EEG:
% $$$ read_EEG_file=@(file_name,channel_names) edfread_nicolet(file_name,channel_names);


%---------------------------------------------------------------------
% 1. load from EDF
%---------------------------------------------------------------------

% if multiple files, then concatenate:
if(iscell(fname))
    N_files=length(fname);
    for n=1:N_files
        fname{n}=strip_file_extension(fname{n});
        [data{n},ch_labels{n},Fs{n}]=read_EEG_file([EEG_DATA_DIR fname{n} '.edf'], ...
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
    [data,ch_labels,Fs]=read_EEG_file([EEG_DATA_DIR fname '.edf'],channel_names);
    
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
% 3. PRE-PROCESS (lowpass filter and resample)
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
isint = @(x) (x == round(x));
if(isint(Fs/Fs_new))
    idec=1:(Fs/Fs_new):size(data,2);
    
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
% 4. SAVE
%---------------------------------------------------------------------
if(SAVE_DATA)
    if(isempty(fname_out))
        if(iscell(fname))
            fname_stub=num2str(cell2mat(regexp(fname{1},'[\d+]','match')));
            mfname=[EEG_DATA_DIR_MATFILES filesep fname_stub '.mat'];
        else
            mfname=[EEG_DATA_DIR_MATFILES filesep fname '.mat'];
        end
    else
        mfname=[EEG_DATA_DIR_MATFILES fname_out '.mat'];
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
% 5. PLOT
%---------------------------------------------------------------------
if(DBplot)
  eeg_plotgui_withannos('signals',eeg_data,'Fs',Fs,'channel_labels',ch_labels,...
                        'bipolar_montage',-1,'epoch_length',60*10);
end


function fname=strip_file_extension(fname)
if(length(fname)>4 && strcmp(fname(end-3:end),'.edf'))
    fname=fname(1:end-4);
end
