%-------------------------------------------------------------------------------
% size_feature: size of feature (e.g. array or 1 point)
%
% Syntax: N=size_feature(feat_name)
%
% Inputs: 
%     feat_name - 
%
% Outputs: 
%     N - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 04-05-2016
%
% last update: Time-stamp: <2016-10-17 09:40:17 (otoolej)>
%-------------------------------------------------------------------------------
function N=size_feature(feat_name)

QUEEN_parameters;

feat_group=strsplit(feat_name,'_');
feat_group=feat_group{1};

if(strcmp(feat_group,'IBI'))
    N=1;
    return;
end


switch feat_name
    %---------------------------------------------------------------------
    % SPECIAL cases
    %---------------------------------------------------------------------
  case 'spectral_edge_frequency'
    N=1;
    
  otherwise
    %---------------------------------------------------------------------
    % feature set is size number of frequency bands
    %---------------------------------------------------------------------
    N=size( feat_params_st.(char(feat_group)).freq_bands,1 );
    
end

