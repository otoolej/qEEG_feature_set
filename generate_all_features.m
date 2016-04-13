%-------------------------------------------------------------------------------
% generate_all_features: generate all features for EEG recording (in .mat)
%
% Syntax: feat_st=generate_all_features(fname,channel_names,feat_set)
%
% Inputs: 
%     fname,channel_names,feat_set - 
%
% Outputs: 
%     feat_st - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 07-04-2016
%
% last update: Time-stamp: <2016-04-13 04:14:59 (otoolej)>
%-------------------------------------------------------------------------------
function feat_st=generate_all_features(fname,channel_names,feat_set)
if(nargin<2 || isempty(channel_names)), channel_names=[]; end
if(nargin<3 || isempty(feat_set)), feat_set=[]; end


quant_feats_parameters;


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
    fprintf('<strong> loading EEG data from file saved on %s </strong>\n',datestr(d.time_now));
    eeg_data=d.eeg_data; Fs=d.Fs; ch_labels=d.ch_labels;
    
end


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


[N_channels,N]=size(eeg_data);

%---------------------------------------------------------------------
% 2. generate features
%---------------------------------------------------------------------
if(isempty(feat_set))
    feat_set=FEATURE_SET_ALL;
end

N_feats=length(feat_set);

% A) iterate over features
for n=1:N_feats
    
    % analysis on a per-channel basis
    if( any(strfind(feat_set{n},'spectral')) || ...
        any(strfind(feat_set{n},'amplitude')) )

        % B) iterate over channels
        feats_channel=[]; x_epochs=[]; 
        for c=1:N_channels
            x_epochs=overlap_epochs(eeg_data(c,:)',Fs,EPOCH_LENGTH,EPOCH_OVERLAP);
            N_epochs=size(x_epochs,1);
            
            % C) iterate over epochs
            feats_epochs=[];
            for e=1:N_epochs
                if(any(strfind(feat_set{n},'spectral')))
                    feats_epochs(e,:)=spectral_features(x_epochs(e,:),Fs,feat_set{n});
                    
                elseif(any(strfind(feat_set{n},'amplitude')))
                    feats_epochs(e,:)=amplitude_features(x_epochs(e,:),Fs,feat_set{n});
                end
            end
            % median over all epochs
            feats_channel(c,:)=nanmedian(feats_epochs);
        end
        % and median over all channels:
        feat_st.(char(feat_set{n}))=nanmedian(feats_channel);

    % or use all channels:
    elseif(strfind(feat_set{n},'connectivity'))

        x_epochs=[]; 
        for c=1:N_channels
            x_epochs(c,:,:)=overlap_epochs(eeg_data(c,:)',Fs,EPOCH_LENGTH,EPOCH_OVERLAP);
        end
        N_epochs=size(x_epochs,2); 
        
        % B) iterate over epochs:
        feats_epochs=[]; x_ep=[];
        for e=1:N_epochs
            x_ep=reshape(x_epochs(:,e,:),size(x_epochs,1),size(x_epochs,3));
            feats_epochs(e,:)=connectivity_features(x_ep,Fs,feat_set{n});
        end
        % median over all epochs
        feat_st.(char(feat_set{n}))=nanmedian(feats_epochs);

    end
end





function [x_epochs]=overlap_epochs(x,Fs,L_window,overlap,window_type)
%---------------------------------------------------------------------
% overlapping epochs in one matrix
%---------------------------------------------------------------------
if(nargin<4 || isempty(overlap)), overlap=50; end
if(nargin<5 || isempty(window_type)), window_type='rect'; end


[L_hop,L_epoch,win_epoch]=get_epoch_window(overlap,L_window,window_type,Fs);

N=length(x);
N_epochs=floor( (N-L_epoch)/L_hop );
if(N_epochs<1) N_epochs=1; end
nw=0:L_epoch-1;


x_epochs=zeros(N_epochs,L_epoch);
for k=1:N_epochs
    nf=mod(nw+(k-1)*L_hop,N);
    
    x_epochs(k,:)=x(nf+1).*win_epoch';
end
