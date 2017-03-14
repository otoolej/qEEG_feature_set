%-------------------------------------------------------------------------------
% set_bbiploar_montage: Convert monopolar (referential) to bi-polar montgage
%
% Syntax: [bi_sigs,bi_labels]=set_bi_montage(sigs,channel_names)
%
% Inputs: 
%     sigs          - EEG data in referential montage
%     channel_names - cell of referential channel names
%                      e.g. {'C3','C4','F3','F4'}
%
% Outputs: 
%     bi_sigs   - EEG data in referential montage
%     bi_labels - cell of bipolar channel names
%                 e.g. {'C3-O1','C4-O2', 'F3-C3', 'F4-C4'}
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%
%     x_ref=data_st.eeg_data_ref;
%     ch_labels_ref=data_st.ch_labels_ref;
%     BI_MONT={{'F4','C4'},{'F3','C3'},{'C4','T4'},{'C3','T3'},{'C4','Cz'}, ...
%              {'Cz','C3'},{'C4','O2'},{'C3','O1'}};
%
%     [x_bi,ch_labels_bi]=set_bi_montage(x_ref,ch_labels_ref,BI_MONT);
%
%     fprintf('bipolar channels: %s\n',strjoin(ch_labels_bi,', '));


% John M. O' Toole, University College Cork
% Started: 05-03-2013
%
% last update: Time-stamp: <2017-03-14 16:18:10 (otoolej)>
%-------------------------------------------------------------------------------
function [bi_sigs,bi_labels]=set_bi_montage(sigs,channel_names,bi_mont)
if(nargin<3), error('requires 3 input arguments.'); end


L=length(bi_mont);  [M,N]=size(sigs);

bi_sigs=zeros(L,N);

for n=1:L
    isig_first=find(~cellfun(@isempty, strfind(upper([channel_names]), ...
                                               upper(bi_mont{n}{1}))));
    isig_second=find(~cellfun(@isempty, strfind(upper([channel_names]), ...
                                                upper(bi_mont{n}{2}))));
    
    bi_sigs(n,:)=sigs(isig_first,:)-sigs(isig_second,:);
end

if(strcmp(channel_names{end},'ECG')==1)
    bi_labels{1}='ECG';
    istart=1;
    bi_sigs(L+1,:)=sigs(L+2,:);
else
    istart=0;
end

for n=1:L
    bi_labels{n+istart}=char( [bi_mont{n}{1} '-' bi_mont{n}{2}] );
end


Le=length(bi_labels);
if(strcmp(channel_names{end},'ECG')==1)
    tmp{1}=bi_labels{1};
    for b=1:L
        tmp{b+1}=bi_labels{Le-b+1};
    end
else
    for b=1:L
        tmp{b}=bi_labels{Le-b+1};
    end
end
bi_labels=tmp;


%bi_labels=bi_labels{end:-1:1};
bi_sigs=bi_sigs(end:-1:1,:);
