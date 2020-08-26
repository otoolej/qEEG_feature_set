%-------------------------------------------------------------------------------
% generate_all_features: generate all features for EEG recording (in .mat)
%
% Syntax: feat_st=generate_all_features(fname,channel_names,feat_set)
%
% Inputs: 
%     fname             - either EEG filename or data structure with EEG 
%                         e.g. data structure: data_st = gen_test_EEGdata(5*60,64,1); 
%     channel_names     - channel labels to process (default, process all)
%                         e.g. {'C3-O1','C4-O2','F3-C3','F4-C4'}
%     feat_set          - cell of features to compute, 
%                         e.g. {'spectral_relative_power','rEEG_SD', 'connectivity_BSI'}
%     return_feat_epoch - return features for each epoch? (true/false; default=false)
%                         (if true then returns 2nd output argument, 'feats_per_epoch')
% 
% Outputs: 
%     feat_st         - structure containing features
%     feats_per_epoch - table of features per epoch; in long format, with columns: 
%                         [time - channel - freq. band - value - feature name]
%                       start_time_sec: start time of the epoch (seconds)
%                       channel: name of channel (string)
%                       freq_band: number of frequency band (string, e.g. FB1)
%                       feature_value: value of the feature (real number, scalar)
%                       feature: name of the feature (string)         
%
% 
% Example:
%       % generate 5 minutes of simulated multichannel EEG data, with 64 Hz sample frequency
%       data_st=gen_test_EEGdata(5*60,64,1);
%
%       % select features to compute:
%       feature_set={'spectral_relative_power','rEEG_SD', 'connectivity_BSI'};
%
%       % generate all features:
%       feat_st=generate_all_features(data_st,[],feature_set);
% 

% John M. O' Toole, University College Cork
% Started: 07-04-2016
%
% last update: Time-stamp: <2020-08-17 18:08:33 (otoolej)>
%-------------------------------------------------------------------------------
function [feat_st, feats_all_epochs_tb] = generate_all_features(fname, channel_names, feat_set, ...
                                                      return_feat_epoch)
if(nargin<2 || isempty(channel_names)), channel_names=[]; end
if(nargin<3 || isempty(feat_set)), feat_set=[]; end
if(nargin<4 || isempty(return_feat_epoch)), return_feat_epoch = false; end



neural_parameters;


%---------------------------------------------------------------------
% 1. load EEG data from .mat file
%---------------------------------------------------------------------
if(isstruct(fname))
    % can input data as structure:
    eeg_data=fname.eeg_data; Fs=fname.Fs; ch_labels=fname.ch_labels;
    
else
    if(length(fname)>4 && strcmp(fname(end-3:end),'.mat'))
        fname=fname(1:end-4);
    end

    % load from .mat file:
    d=load([EEG_DATA_DIR_MATFILES fname '.mat']);
    fprintf(col_str(' loading EEG data from file saved on %s\n',1),datestr(d.time_now));
    eeg_data=d.eeg_data; Fs=d.Fs; ch_labels=d.ch_labels;

end
feats_per_epochs=[];


% select channels:
if(~isempty(channel_names))
    ikeep=[];
    for n=1:length(channel_names)
        it=find( strcmp(ch_labels,channel_names{n}) );
        if(~isempty(it))
            ikeep=[ikeep it];
        end
    end
    eeg_data=eeg_data(ikeep,:);
    ch_labels=ch_labels(ikeep);
end

% or remove empty channels:
irem=[];
for n=1:length(ch_labels)
    if( all(isnan(eeg_data(n,:))) )
        irem=[irem n];
    end
end
if(~isempty(irem))
% $$$     fprintf('removing channels: %s\n',ch_labels(irem));
    eeg_data(irem,:)=[]; ch_labels(irem)=[];
end


[N_channels,N]=size(eeg_data);

%---------------------------------------------------------------------
% 2. generate features
%---------------------------------------------------------------------
if(isempty(feat_set))
    feat_set=FEATURE_SET_ALL;
end

N_feats=length(feat_set);
feats_all_epochs_tb = [];

% A) iterate over features
for n=1:N_feats

    L_feature=size_feature(feat_set{n});
    feat_group=strsplit(feat_set{n},'_');
    feat_group=feat_group{1};

    %---------------------------------------------------------------------
    % SPECTRAL and AMPLITUDE
    % (analysis on a per-channel basis and divide each channel into epochs)
    %---------------------------------------------------------------------
    if( any(strcmp({'amplitude','spectral','rEEG','FD'},feat_group)) )

        % B) iterate over channels
        feats_channel=[]; x_epochs=[]; feats_tbl = [];
        for c=1:N_channels
            [x_epochs, epoch_start_times] = overlap_epochs(eeg_data(c,:)',Fs,EPOCH_LENGTH,EPOCH_OVERLAP);
            N_epochs=size(x_epochs,1);
            
            % C) iterate over epochs
            feats_epochs=NaN(N_epochs,L_feature);
            for e=1:N_epochs
                L_nans=length(find(isnan(x_epochs(e,:))));
                
                if(100*(L_nans/length(x_epochs(e,:))) < EPOCH_IGNORE_PRC_NANS)
                    if(strcmp(feat_group,'spectral'))
                        feats_epochs(e,:)=spectral_features(x_epochs(e,:),Fs, ...
                                                            feat_set{n});
                        
                    elseif(strcmp(feat_group,'FD'))
                        feats_epochs(e,:)=fd_features(x_epochs(e,:),Fs);
                        
                    elseif(strcmp(feat_group,'amplitude'))
                        feats_epochs(e,:)=amplitude_features(x_epochs(e,:),Fs, ...
                                                             feat_set{n});
                        
                    elseif(strcmp(feat_group,'rEEG'))
                        feats_epochs(e,:)=rEEG(x_epochs(e,:),Fs,feat_set{n});
                        
                    end
                end
            end
            % if want to return feature estimated over all epochs:
            if(return_feat_epoch)
        	% feats_per_epochs{n}(c,:,:)=feats_epochs;
                
                % create table with features and start time of epoch:
                fb_names = arrayfun(@(x) ['FB' num2str(x)], 1:size(feats_epochs, 2), 'un', false);
                tb = array2table([epoch_start_times' feats_epochs], ...
                                 'VariableNames', ['start_time_sec', fb_names]);
                % add channel:
                tb.channel(:) = string(ch_labels{c});
                
                % convert from wide to long format for frequency bands:
                tb = stack(tb, fb_names, 'newDataVariableName', {'feature_value'}, ...
                           'IndexVariableName', {'freq_band'});
                
                feats_tbl = [feats_tbl; tb];
            end
            
            % median over all epochs
            feats_channel(c,:)=nanmedian(feats_epochs, 1);
        end
        % and median over all channels:
        feat_st.(char(feat_set{n}))=nanmedian(feats_channel, 1);
        

        if(return_feat_epoch)
            % add feature name and combine:
            feats_tbl.feature(:) = string(feat_set{n});
            feats_all_epochs_tb = [feats_all_epochs_tb; feats_tbl];
        end


        %---------------------------------------------------------------------
        % CONNECTIVITY FEATURES
        % (use over all channels but also divide into epochs)
        %---------------------------------------------------------------------
    elseif(strfind(feat_set{n},'connectivity'))

        x_epochs=[]; 
        for c=1:N_channels
            if(c == N_channels)
                [x_epochs(c,:,:), epoch_start_times] = ...
                    overlap_epochs(eeg_data(c,:)',Fs,EPOCH_LENGTH,EPOCH_OVERLAP);
            else
                x_epochs(c,:,:) = overlap_epochs(eeg_data(c,:)',Fs,EPOCH_LENGTH,EPOCH_OVERLAP);
            end
        end
        N_epochs=size(x_epochs,2); 
        
        % B) iterate over epochs:
        feats_epochs=NaN(N_epochs,L_feature);
        x_ep=[];
        
        for e=1:N_epochs
            x_ep=reshape(x_epochs(:,e,:),size(x_epochs,1),size(x_epochs,3));
            
            L_nans=length(find(isnan(x_ep(:))));
            if(100*(L_nans/length(x_ep(:))) < EPOCH_IGNORE_PRC_NANS)
                
                feats_epochs(e,:)=connectivity_features(x_ep,Fs,feat_set{n},[], ...
                                                        ch_labels);
            end
            
        end
        % median over all epochs
        feat_st.(char(feat_set{n}))=nanmedian(feats_epochs, 1);
        
        % if want to return feature estimated over all epochs:
        if(return_feat_epoch)
            % create table with features and start time of epoch:
            fb_names = arrayfun(@(x) ['FB' num2str(x)], 1:size(feats_epochs, 2), 'un', false);
            tb = array2table([epoch_start_times' feats_epochs], ...
                             'VariableNames', ['start_time_sec', fb_names]);
            % add channel:
            tb.channel(:) = NaN;
            
            % convert from wide to long format for frequency bands:
            tb = stack(tb, fb_names, 'newDataVariableName', {'feature_value'}, ...
                       'IndexVariableName', {'freq_band'});
            
            % add feature name and combine:
            tb.feature(:) = string(feat_set{n});
            feats_all_epochs_tb = [feats_all_epochs_tb; tb];
        end
        
        

        %---------------------------------------------------------------------
        % inter-burst interval features
        % (use entire recording but channel-by-channel)
        %---------------------------------------------------------------------
    elseif(strfind(feat_set{n},'IBI_'))
        
        % B) iterate over channels
        feats_channel=NaN(N_channels,L_feature); 
        for c=1:N_channels
            feats_channel(c,:)=IBI_features(eeg_data(c,:)',Fs,feat_set{n});
        end
        % and median over all channels:
        feat_st.(char(feat_set{n}))=nanmedian(feats_channel,1);


    end
end




function [x_epochs, start_times] = overlap_epochs(x, Fs, L_window, overlap, window_type)
%---------------------------------------------------------------------
% overlapping epochs in one matrix
%---------------------------------------------------------------------
if(nargin<4 || isempty(overlap)), overlap=50; end
if(nargin<5 || isempty(window_type)), window_type='rect'; end

x = x(:).';

[L_hop, L_epoch, win_epoch] = gen_epoch_window(overlap, L_window, window_type, Fs);

N = length(x);
N_epochs = ceil( (N - (L_epoch - L_hop)) / L_hop );
if(N_epochs < 1) 
    N_epochs = 1; 
    fprintf('| WARNING: signal length is less than segment length (L_epoch - L_hop).\n');
    fprintf('| Adjust ''EPOCH_LENGTH'' or ''EPOCH_OVERLAP'' in ');
    fprintf('''neural_parameters.m'' file.\n');
end
nw = 0:(L_epoch - 1);
ix = 0:(N - 1);


x_epochs = NaN(N_epochs, L_epoch);
start_times = NaN(1, N_epochs);
for k = 1:N_epochs
    nf = nw + (k - 1) * L_hop;
    % zero-pad if outside x:
    nf = nf(ismember(nf, ix)) + 1;
    i_nf = 1:length(nf);
    start_times(k) = (min(nf) - 1) / Fs;

    x_epochs(k, i_nf) = x(nf) .* win_epoch(i_nf).';
end
