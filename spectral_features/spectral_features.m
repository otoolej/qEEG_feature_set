%-------------------------------------------------------------------------------
% gen_spectral_power: spectral power for all frequency bands
%
% Syntax: featx=spectral_features(x,Fs,feat_name,params_st)
%
% Inputs: 
%     x          - epoch of EEG data (size 1 x N)
%     Fs         - sampling frequency (in Hz)
%     feat_name  - feature type, defaults to 'spectral_power';
%                  see full list of 'spectral_' features in all_features_list.m
%     params_st  - parameters (as structure); 
%                  see neural_parameters.m for examples
%
% Outputs: 
%     featx  - feature at each frequency band 
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%
%     featx=spectral_features(x,Fs,'spectral_relative_power');
%     
%

% John M. O' Toole, University College Cork
% Started: 07-04-2016
%
% last update: Time-stamp: <2020-08-17 16:53:32 (otoolej)>
%-------------------------------------------------------------------------------
function featx=spectral_features(x,Fs,feat_name,params_st)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='spectral_power'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end

DBplot=0;

if(isempty(params_st))
    neural_parameters;
    if(strfind(feat_name,'spectral'))
        params_st=feat_params_st.spectral;
    else
        params_st=feat_params_st.(char(feat_name));
    end
end

freq_bands=params_st.freq_bands;
total_freq_bands=params_st.total_freq_bands;


if(length(x) < (params_st.L_window * Fs))
    warning('SPECTRAL features: signal length < window length; set shorter L_window');
    warning('SPECTRAL features: not calculating spectral features.');
    featx = NaN;
    return;
end


switch feat_name
  case {'spectral_power','spectral_relative_power'}
    %---------------------------------------------------------------------
    % use periodogram to estimate spectral power 
    %---------------------------------------------------------------------
    params_st.method='periodogram';
    [pxx,itotal_bandpass,f_scale,N,fp]=gen_spectrum(x,Fs,params_st,1);
    pxx=pxx.*Fs;

    Nh=length(pxx);
    
    if(DBplot)
        figure(1); clf; hold all;
        plot(fp,20*log10(pxx));
    end
    
    if(strcmp(feat_name,'spectral_relative_power'))
        pxx_total=sum( pxx(itotal_bandpass) )/N;
    else
        pxx_total=1;
    end
    
    spec_pow=NaN(1,size(freq_bands,1));

    for p=1:size(freq_bands,1)
        if(p==1)
            istart=ceil(freq_bands(p,1)*f_scale);
        else
            istart=ibandpass(end)-1;
        end
        ibandpass=istart:floor(freq_bands(p,2)*f_scale);        
        ibandpass=ibandpass+1;
        ibandpass(ibandpass<1)=1; ibandpass(ibandpass>Nh)=Nh;    
        
        spec_pow(p)=sum( pxx(ibandpass) )/(N*pxx_total);            
        
        if(DBplot)
            line([fp(ibandpass(1)) fp(ibandpass(1))],ylim,'color','k');        
            line([fp(ibandpass(end)) fp(ibandpass(end))],ylim,'color','k');
        end
    end
    featx=spec_pow;
    
        
  case 'spectral_flatness'
    %---------------------------------------------------------------------
    % spectral flatness (Wiener entropy)
    %---------------------------------------------------------------------
    [pxx,~,f_scale]=gen_spectrum(x,Fs,params_st);


    % for each frequency band:
    N_freq_bands=size(freq_bands,1);
    featx=NaN(1,N_freq_bands);
    N=length(pxx);    
    for p=1:N_freq_bands
        if(p==1)
            istart=ceil(freq_bands(p,1)*f_scale);
        else
            istart=ibandpass(end)-1;
        end
        ibandpass=istart:floor(freq_bands(p,2)*f_scale);        
        ibandpass=ibandpass+1;
        ibandpass(ibandpass<1)=1; ibandpass(ibandpass>N)=N;    
        
        % geometric mean / arthimetric mean:
        featx(p)=exp(nanmean(log(pxx(ibandpass)+eps)))/nanmean(pxx(ibandpass));
    end
    

  case 'spectral_entropy'
    %---------------------------------------------------------------------
    % spectral entropy (= Shannon entropy on normalised PSD)
    %---------------------------------------------------------------------
    [pxx,~,f_scale]=gen_spectrum(x,Fs,params_st);


    % for each frequency band:
    N_freq_bands=size(freq_bands,1);
    featx=NaN(1,N_freq_bands);
    N=length(pxx);    
    for p=1:N_freq_bands
        if(p==1)
            istart=ceil(freq_bands(p,1)*f_scale);
        else
            istart=ibandpass(end)-1;
        end
        ibandpass=istart:floor(freq_bands(p,2)*f_scale);        
        ibandpass=ibandpass+1;
        ibandpass(ibandpass<1)=1; ibandpass(ibandpass>N)=N;    

        pr=pxx(ibandpass)./sum(pxx(ibandpass));
    
        featx(p)=-sum( pr.*log(pr+eps) )/log(length(pr));
    end    

    
  case 'spectral_diff'
    %---------------------------------------------------------------------
    % spectral difference using the spectrogram
    %---------------------------------------------------------------------
    % a) generate spectrogram
    [S_stft,~,f_scale]=gen_STFT(x,params_st.L_window,params_st.window_type, ...
                                params_st.overlap,Fs);
    
    [N_epochs,M]=size(S_stft);

    N_freq_bands=size(freq_bands,1);
    featx=NaN(1,N_freq_bands);
    for p=1:N_freq_bands
        if(p==1)
            istart=ceil(freq_bands(p,1)*f_scale);
        else
            istart=ibandpass(end)-1;
        end
        ibandpass=istart:floor(freq_bands(p,2)*f_scale);        
        ibandpass=ibandpass+1;
        ibandpass(ibandpass<1)=1; ibandpass(ibandpass>M)=M;    

        S_stft_band=S_stft(:,ibandpass)./max(max(S_stft(:,ibandpass)));

        spec_diff=zeros(1,N_epochs - 1);
        for n=1:N_epochs-1
            v1=S_stft_band(n,:);  v2=S_stft_band(n+1,:);
            spec_diff(n)=mean( abs(v1-v2).^2 );
        end

        % how to summarise this?? 
        featx(p)=nanmedian(spec_diff);

    % or return continuous function (need to interpolate to write sample frequency):
% $$$     sd_int=interp1( (1:N_epochs)./(N_epochs/N), spec_diff, 1:N, 'linear' );
% $$$     featx(p,:)=sd_int;
    end
    
  case 'spectral_edge_frequency'
    %---------------------------------------------------------------------
    % spectral edge frequency
    %---------------------------------------------------------------------
    [pxx,itotal_bandpass,~,~,fp]=gen_spectrum(x,Fs,params_st);
    
    % only within this frequency band:
    pxx(setdiff(1:length(pxx),itotal_bandpass))=0;

    pxx=pxx./sum(pxx);

    % compute cumulative density
    pxx_cum=cumsum(pxx);

    % spectral edge frequency corresponds to the frequency at (nearest to)
    % the point on the freq axis where pyy_cum = 0.05
    [vv, idx]=min(abs(pxx_cum-params_st.SEF));
    featx=fp(idx);

    
    if(DBplot)
        figure(2); clf; hold all;
        hx(1)=subplot(211); hold all;
        plot(fp,pxx_cum);
        line([1 1].*featx,ylim,'color','r');
        line(xlim,[1 1].*params_st.SEF,'color','r');        
        hx(1)=subplot(212); hold all;
        plot(fp,10*log10(pxx));
    end
    
  otherwise
    fprintf('unknown feature ''%s''; check spelling\n',feat_name);
    featx=NaN;
end








end

