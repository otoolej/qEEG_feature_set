%-------------------------------------------------------------------------------
% amplitude_features: amplitude features
%
% Syntax: featx=amplitude_features(x,Fs,feat_name,params_st)
%
% Inputs: 
%     x,Fs,feat_name,params_st - 
%
% Outputs: 
%     featx - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 12-04-2016
%
% last update: Time-stamp: <2016-11-03 15:51:59 (otoolej)>
%-------------------------------------------------------------------------------
function featx=amplitude_features(x,Fs,feat_name,params_st)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='envelope'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end

DBplot=0;

% so far no parameters but maybe later:
if(isempty(params_st))
    neural_parameters;
    if(strfind(feat_name,'amplitude'))
        params_st=feat_params_st.amplitude;
    else
        params_st=feat_params_st.(char(feat_name));
    end
end
freq_bands=params_st.freq_bands;

N_freq_bands=size(freq_bands,1);
if(isempty(freq_bands))
    N_freq_bands=1;
end

x_orig=x;

for n=1:N_freq_bands
    
    if(~isempty(freq_bands))
        x=filter_butterworth_withnans(x_orig,Fs,freq_bands(n,2),freq_bands(n,1),5, ...
                                      params_st.FILTER_REPLACE_ARTEFACTS);
    end
    

    switch feat_name
      case 'amplitude_total_power'
        %---------------------------------------------------------------------
        % power of signal
        %---------------------------------------------------------------------
        featx(n)=nanmean( abs(x).^2 );
        
      case {'amplitude_env_mean','amplitude_env','amplitude_env_SD'}
        %---------------------------------------------------------------------
        % mean or SD of envelope
        %---------------------------------------------------------------------
        x(isnan(x))=[];
        env=abs( hilbert(x) ).^2;
        
        if(strcmp(feat_name,'amplitude_env_mean') || strcmp(feat_name,'amplitude_env'))
            featx(n)=nanmean(env);
        elseif(strcmp(feat_name,'amplitude_env_SD'))
            featx(n)=std(env);
        end
        
      case 'amplitude_SD'
        %---------------------------------------------------------------------
        % standard deviation of amplitude
        %---------------------------------------------------------------------
        featx(n)=nanstd(x);
        
      case 'amplitude_skew'
        %---------------------------------------------------------------------
        % skew of amplitude
        %---------------------------------------------------------------------
        featx(n)=abs(skewness(x));
        
      case 'amplitude_kurtosis'
        %---------------------------------------------------------------------
        % kurtosis of amplitude
        %---------------------------------------------------------------------
        featx(n)=kurtosis(x);
        
      otherwise
        error(['unknown feature: ' feat_name]);
    end
    
    
    
end


