%-------------------------------------------------------------------------------
% gen_spectral_power: spectral power for all frequency bands
%
% Syntax: featx=spectral_features(x,Fs)
%
% Inputs: 
%     x,Fs - 
%
% Outputs: 
%     featx - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 07-04-2016
%
% last update: Time-stamp: <2016-04-14 18:16:02 (otoolej)>
%-------------------------------------------------------------------------------
function featx=spectral_features(x,Fs,feat_name,params_st)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='spectral_power'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end

DBplot=0;

if(isempty(params_st))
    quant_feats_parameters;
    if(strfind(feat_name,'spectral'))
        params_st=feat_params_st.spec;
    else
        params_st=feat_params_st.(char(feat_name));
    end
end

freq_bands=params_st.freq_bands;
total_freq_bands=params_st.total_freq_bands;

switch feat_name
  case {'spectral_power','relative_spectral_power'}

    % 2 different options for generating spectral power measures:
    switch params_st.method
      case 'PSD'
        %---------------------------------------------------------------------
        % Welch's PSD estimate
        %---------------------------------------------------------------------
        [pxx,itotal_bandpass,f_scale,fp]=psd_Welch(x,params_st.L_window,params_st.window_type, ...
                                                   params_st.overlap,freq_bands,total_freq_bands,Fs);
        
        if(DBplot)
            figure(1); clf; hold all;
            plot(fp,20*log10(pxx));
        end
        
        if(strcmp(feat_name,'relative_spectral_power'))
            pxx_total=sum( pxx(itotal_bandpass) );
        else
            pxx_total=1;
        end
        
        spec_pow=NaN(1,size(freq_bands,1));
        N=length(pxx);

        for p=1:size(freq_bands,1)
            ibandpass=ceil(freq_bands(p,1)*f_scale):floor(freq_bands(p,2)*f_scale);        
            ibandpass=ibandpass+1;
            ibandpass(ibandpass<1)=1; ibandpass(ibandpass>N)=N;    
            
            
            spec_pow(p)=sum( pxx(ibandpass) )/pxx_total;
            
            if(DBplot)
                line([fp(ibandpass(1)) fp(ibandpass(1))],ylim,'color','k');        
                line([fp(ibandpass(end)) fp(ibandpass(end))],ylim,'color','k');
            end
        end
        featx=spec_pow;
        

      case {'spectrogram','med-spectrogram'}
        %---------------------------------------------------------------------
        % estimate spectral power on short-time FFT and then average
        %---------------------------------------------------------------------
        [S_stft,Nfreq,f_scale]=gen_STFT(x,params_st.L_window,params_st.window_type,...
                                        params_st.overlap,Fs);
        
        [N_epochs,M]=size(S_stft);

        if(strcmp(feat_name,'relative_spectral_power'))
            itotal_bandpass=ceil(total_freq_bands(1)*f_scale):floor(total_freq_bands(2)*f_scale);
            itotal_bandpass=itotal_bandpass+1;
            itotal_bandpass(itotal_bandpass<1)=1; itotal_bandpass(itotal_bandpass>M)=M;  
        end
        
        spec_pow=NaN(N_epochs,size(freq_bands,1));
        
        for p=1:size(freq_bands,1)
            ibandpass=ceil(freq_bands(p,1)*f_scale):floor(freq_bands(p,2)*f_scale);        
            ibandpass=ibandpass+1;
            ibandpass(ibandpass<1)=1; ibandpass(ibandpass>M)=M;    

            
            for k=1:N_epochs
                pxx=S_stft(k,:);
                
                if(strcmp(feat_name,'relative_spectral_power'))
                    spec_pow(k,p)=sum( pxx(ibandpass) )/sum( pxx(itotal_bandpass) );
                else
                    spec_pow(k,p)=sum( pxx(ibandpass) );
                end
            end
        end
        
        if(strcmp(params_st.method,'med-spectrogram'))
            featx=nanmedian(spec_pow);
        else
            featx=nanmean(spec_pow);
        end

    end
    
    
  case 'spectral_flatness'
    %---------------------------------------------------------------------
    % spectral flatness (Wiener entropy)
    %---------------------------------------------------------------------
    [pxx,itotal_bandpass]=psd_Welch(x,params_st.L_window,params_st.window_type, ...
                                    params_st.overlap, freq_bands,total_freq_bands,Fs);
    pxx=pxx(itotal_bandpass);
    
    N=length(pxx);
    
    % geometric mean / arthimetric mean:
    featx=exp(nanmean(log(pxx+eps)))/nanmean(pxx);
    
% $$$     featx=(prod(pxx)^(1/N))/nanmean(pxx);      
    

  case 'spectral_entropy'
    %---------------------------------------------------------------------
    % spectral entropy (= Shannon entropy on normalised PSD)
    %---------------------------------------------------------------------
    [pxx,itotal_bandpass]=psd_Welch(x,params_st.L_window,params_st.window_type, ...
                                    params_st.overlap, freq_bands,total_freq_bands,Fs);
    pxx=pxx(itotal_bandpass);

    p=pxx./sum(pxx);
    
    featx=-sum( p.*log(p+eps) )/log(length(pxx));
    
    
  case 'spectral_diff'
    %---------------------------------------------------------------------
    % spectral difference using the spectrogram
    %---------------------------------------------------------------------
    N=length(x); 
    
    % a) generate spectrogram
    [S_stft,Nfreq,f_scale]=gen_STFT(x,params_st.L_window,params_st.window_type,...
                                    params_st.overlap,Fs);
    
    [N_epochs,M]=size(S_stft);
        
    itotal_bandpass=ceil(total_freq_bands(1)*f_scale):floor(total_freq_bands(2)*f_scale);
    itotal_bandpass=itotal_bandpass+1;
    itotal_bandpass(itotal_bandpass<1)=1; itotal_bandpass(itotal_bandpass>M)=M;    

    % normalize to whole spectrogram:
    S_stft=S_stft(:,itotal_bandpass)./max(max(S_stft(:,itotal_bandpass)));

    spec_diff=zeros(1,N_epochs);

    for n=1:N_epochs-1
        v1=S_stft(n,:);  v2=S_stft(n+1,:);
        spec_diff(n)=mean( abs(v1-v2).^2 );
    end

    % how to summarise this?? 
    featx=nanmedian(spec_diff);

    % or return continuous function (need to interpolate to write sample frequency):
% $$$     sd_int=interp1( (1:N_epochs)./(N_epochs/N), spec_diff, 1:N, 'linear' );
% $$$     featx=sd_int;

    
  case 'spectral_edge_frequency'
    %---------------------------------------------------------------------
    % spectral edge frequency
    %---------------------------------------------------------------------
    
    [pxx,itotal_bandpass,~,fp]=psd_Welch(x,params_st.L_window,params_st.window_type, ...
                                    params_st.overlap, freq_bands, total_freq_bands,Fs);
    
    % only within this frequency band:
    pxx(setdiff(1:length(pxx),itotal_bandpass))=0;

    pxx=pxx./sum(pxx);

    % compute cumulative density
    pxx_cum=cumsum(pxx);

    
    % spectral edge frequency corresponds to the frequency at (nearest to)
    % the point on the freq axis where pyy_cum = 0.05
    [vv, idx]=min(abs(pxx_cum-params_st.SEF));
    featx=fp(idx);

    
end




function [pxx,itotal_bandpass,f_scale,fp]=psd_Welch(x,L_window,window_type,overlap, ...
                                                  freq_bands,total_freq_bands,Fs)
%---------------------------------------------------------------------
% generate PSD using the Welch method
%---------------------------------------------------------------------

% a) PSD estimate (Welch's periodogram):
win_length=make_odd(L_window*Fs);
overlap=ceil(win_length*(1-overlap/100));

[pxx,fp]=pwelch(x,win_length,overlap,[],Fs);

% b) limit to frequency band of interest:
N=length(pxx);
Nfreq=2*(N-1);
f_scale=(Nfreq/Fs);
itotal_bandpass=ceil(total_freq_bands(1)*f_scale):floor(total_freq_bands(2)*f_scale);
itotal_bandpass=itotal_bandpass+1;
itotal_bandpass(itotal_bandpass<1)=1;  itotal_bandpass(itotal_bandpass>N)=N;    




function [S_stft,Nfreq,f_scale]=gen_STFT(x,L_window,window_type,overlap,Fs)
%---------------------------------------------------------------------
% Short-time Fourier transform (magnitude only, i.e. spectrogram):
%---------------------------------------------------------------------
[L_hop,L_epoch,win_epoch]=get_epoch_window(overlap,L_window,window_type,Fs);

N=length(x);
N_epochs=floor( (N-L_epoch)/L_hop );
if(N_epochs<1) N_epochs=1; end
nw=0:L_epoch-1;
Nfreq=2^nextpow2(L_epoch);

%---------------------------------------------------------------------
% generate short-time FT on all data:
%---------------------------------------------------------------------
K_stft=zeros(N_epochs,L_epoch);
for k=1:N_epochs
    nf=mod(nw+(k-1)*L_hop,N);
    
    K_stft(k,:)=x(nf+1).*win_epoch';
end
S_stft=abs(fft(K_stft.',Nfreq)).^2;
S_stft=S_stft';


% return only positive frequencies:
S_stft=S_stft(:,1:ceil(Nfreq/2)+1);

f_scale=(Nfreq/Fs);
