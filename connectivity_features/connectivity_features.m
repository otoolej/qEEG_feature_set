%-------------------------------------------------------------------------------
% connectivity_features: connectivity features
%
% Syntax: featx=connectivity_features(x,Fs,feat_name,params_st)
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
% Started: 13-04-2016
%
% last update: Time-stamp: <2016-04-14 18:17:31 (otoolej)>
%-------------------------------------------------------------------------------
function featx=connectivity_features(x,Fs,feat_name,params_st)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='envelope'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end

DBplot=0;

% so far no parameters but maybe later:
if(isempty(params_st))
    quant_feats_parameters;
    if(strfind(feat_name,'connectivity'))
        params_st=feat_params_st.connectivity;
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

N=size(x,2);
Nfreq_h=ceil(N/2);
f_scale=(N/Fs);
itotal_bandpass=ceil(freq_bands(1,1)*f_scale):floor(freq_bands(end,2)*f_scale);
itotal_bandpass=itotal_bandpass+1;
itotal_bandpass(itotal_bandpass<1)=1;  itotal_bandpass(itotal_bandpass>Nfreq_h)=Nfreq_h;    


for n=1:N_freq_bands
    
    ibandpass=ceil(freq_bands(n,1)*f_scale):floor(freq_bands(n,2)*f_scale);        
    ibandpass=ibandpass+1;
    ibandpass(ibandpass<1)=1; ibandpass(ibandpass>Nfreq_h)=Nfreq_h;    
        

    switch feat_name
      case 'connectivity_BSI'
        %---------------------------------------------------------------------
        % brain sysmetry index (Van Putten, 2007)
        %---------------------------------------------------------------------
        
        % assuming just two pairs corresponding to left/right:
        X_left=abs(fft(x(1,:))).^2; 
        X_right=abs(fft(x(2,:))).^2;

        X_left=X_left(ibandpass);
        X_right=X_right(ibandpass);        

        featx(n)=nanmean(abs( (X_left - X_right)./(X_left + X_right) ));
        
        
      case 'connectivity_corr'
        %---------------------------------------------------------------------
        % cross-correlation (Spearmans)
        %---------------------------------------------------------------------
        x_filt(1,:)=filt_butterworth(x(1,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        x_filt(2,:)=filt_butterworth(x(2,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        
        featx(n)=abs(corr(x_filt(1,:)',x_filt(2,:)','type','spearman'));

      case 'connectivity_lag_corr'
        %---------------------------------------------------------------------
        % lag of cross-correlation 
        %---------------------------------------------------------------------
        x_filt(1,:)=filt_butterworth(x(1,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        x_filt(2,:)=filt_butterworth(x(2,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        
        [cc,lag]=xcorr(x_filt(1,:),x_filt(2,:),'biased');

        [~,imax]=max(abs(cc));
        time_max_lag=lag(imax)/Fs;

% $$$         figure(1); clf; hold all;
% $$$         plot(lag,cc);
        
        featx(n)=time_max_lag;
        
        
      case 'connectivity_asynchrony'
        %---------------------------------------------------------------------
        % inter-hemispheric asynchrony (Korotchikova et al. 2011)
        %---------------------------------------------------------------------
        
        % A) BSI:
        X_left=abs(fft(x(1,:))).^2; X_right=abs(fft(x(2,:))).^2;
        X_left=X_left(ibandpass);   X_right=X_right(ibandpass);        

        bsi=1-nanmean(abs( (X_left - X_right)./(X_left + X_right) ));

        % B) correlation coefficient:
        x_filt(1,:)=filt_butterworth(x(1,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        x_filt(2,:)=filt_butterworth(x(2,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        
        cc=abs(corr(x_filt(1,:)',x_filt(2,:)'));

% $$$         dispVars(cc,bsi);
        featx(n)=(cc-bsi);
    end
    
    
    

end
