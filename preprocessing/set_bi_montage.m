%-------------------------------------------------------------------------------
% set_bbiploar_montage: Convert monopolar (referential) to bi-polar montgage
%
% Syntax: [bi_sigs,bi_labels]=set_bi_montage(sigs,channel_names)
%
% Inputs: 
%     sigs,channel_names - 
%
% Outputs: 
%     [bi_sigs,bi_labels] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 05-03-2013
%
% last update: Time-stamp: <2016-04-05 12:58:36 (otoolej)>
%-------------------------------------------------------------------------------
function [bi_sigs,bi_labels]=set_bi_montage(sigs,channel_names,bi_mont)
if(nargin<3), error('requires 3 input arguments.'); end


L=length(bi_mont);  [M,N]=size(sigs);

bi_sigs=zeros(L,N);

for n=1:L
    isig_first=find(~cellfun(@isempty, strfind([channel_names],bi_mont{n}{1})));
    isig_second=find(~cellfun(@isempty, strfind([channel_names],bi_mont{n}{2})));
    
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
