%-------------------------------------------------------------------------------
% IBI_features: burst/inter-burst interval features (assuming preterm EEG, <32 weeks GA)
%
% Syntax: featx=IBI_features(x,Fs,feat_name,params_st)
%
% Inputs: 
%     x          - epoch of EEG data (size 1 x N)
%     Fs         - sampling frequency (in Hz)
%     feat_name  - feature type, defaults to 'IBI_length_max';
%                  see full list of 'IBI_' features in all_features_list.m
%     params_st  - parameters (as structure); 
%                  see neural_parameters.m for examples
%
% Outputs: 
%     featx  - IBI features
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%
%     featx=IBI_features(x,Fs,'IBI_burst_number');
%     
%
% [1] JM O’ Toole, GB Boylan, RO Lloyd, RM Goulding, S Vanhatalo, NJ Stevenson
% (2017). Detecting Bursts in the EEG of Very and Extremely Premature Infants using a
% Multi-Feature Approach. under review, 2017.


% John M. O' Toole, University College Cork
% Started: 26-04-2016
%
% last update: Time-stamp: <2020-04-22 18:44:22 (otoolej)>
%-------------------------------------------------------------------------------
function featx=IBI_features(x,Fs,feat_name,params_st)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='IBI_length_max'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end

DBplot=0;


bdetect_path=exist('eeg_burst_detector.m','file');
if(bdetect_path~=2)
    fprintf('\n** ------------ **\n');    
    fprintf('Burst detector not included in path (eeg_burst_detector.m)\n\n');
    fprintf('If installed, ensure eeg_burst_detector.m and associated files \n');
    fprintf('are included in the Matlab path.\n');
    fprintf('To do this, run the file load_curdir.m in the burst_detector folder.\n');
    fprintf('(for more on search paths see: ');
    fprintf(['<a href=http://uk.mathworks.com/help/matlab/matlab_env/what-is-the-' ...
             'matlab-search-path.html>Matlab search path</a>)\n\n']);
    fprintf('If the burst detector is not installed, download from:\n');
    fprintf(['<a href=http://otoolej.github.io/code/burst_detector/>burst detector ' ...
             'source code </a>\n']);
    fprintf('** ------------ **\n');    
    
    featx=NaN;
    return;
end



[burst_anno,t_stat]=eeg_burst_detector(x,Fs);


if(DBplot)
    figure(1); clf; hold all;
    hx(1)=subplot(211); plot(x);
    hx(2)=subplot(212); plot(burst_anno); ylim([-0.2 1.2])
    linkaxes(hx,'x');
end


%---------------------------------------------------------------------
% remove NaNs
% 
% NOTE -- thisassumes that the same state (either burst or IBI) continues
% after the sequence of NaNs
%---------------------------------------------------------------------
burst_anno = burst_anno(~isnan(burst_anno));    


switch feat_name
  case 'IBI_length_max'
    %---------------------------------------------------------------------
    % max. (95th percentile) inter-burst interval
    %---------------------------------------------------------------------
    featx=estimate_IBI_lengths(burst_anno,95,Fs);    

  case 'IBI_length_median'
    %---------------------------------------------------------------------
    % median inter-burst interval
    %---------------------------------------------------------------------
    featx=estimate_IBI_lengths(burst_anno,50,Fs);    
    
    
  case 'IBI_burst_prc'
    %---------------------------------------------------------------------
    % percentage of bursts
    %---------------------------------------------------------------------
    featx=( length(find([burst_anno]==1))/length([burst_anno]) )*100;
    
    
  case 'IBI_burst_number'
    %---------------------------------------------------------------------
    % number of bursts
    %---------------------------------------------------------------------
    lens_anno=len_zeros(burst_anno,1);    
    featx=length(lens_anno);

  otherwise
    fprintf('unknown feature ''%s''; check spelling\n',feat_name);
    featx=NaN;
    
end



function pc_anno = estimate_IBI_lengths(anno, percentiles_all, Fs)
%---------------------------------------------------------------------
% estimate max./median IBI length
%---------------------------------------------------------------------
min_ibi_interval = Fs / 4;

lens_anno = len_zeros(anno, 0);

ishort = find(lens_anno < min_ibi_interval);
if(~isempty(ishort))
    lens_anno(ishort) = []; 
end
if(~isempty(lens_anno))
    pc_anno = prctile(lens_anno,percentiles_all);
else
    pc_anno = zeros(1, length(percentiles_all));
end
pc_anno = pc_anno ./ Fs;




function [lens,istart,iend]=len_zeros(x,const)
%---------------------------------------------------------------------
% length of continuous runs of zeros
%---------------------------------------------------------------------
if(nargin<2 || isempty(const)), const=0; end

DBplot=0;

x=x(:).';

if( ~all(ismember(sort(unique(x(~isnan(x)))),[0 1])) || ...
    ~ismember(const,[0 1]) )
    warning('must be binary signal');
    return;
end
if(const==1)
    y=x;
    y(~isnan(x))=~x(~isnan(x));
else
    y=x;
end

% find run of zeros:
iedge=diff([0 y==0 0]);
istart=find(iedge==1);
iend=find(iedge==-1)-1;
lens=[iend-istart];

