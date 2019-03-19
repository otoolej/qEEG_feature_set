%-------------------------------------------------------------------------------
% gen_spectrum: specturm using either:
%                             1. periodogram
%                             2. Welch's PSD
%                             3. robust Welch PSD
%
%
% Syntax: [pxx,itotal_bandpass,f_scale,Nfreq,fp]=gen_spectrum(x,Fs,param_st,SCALE_PSD)
%
% Inputs: 
%     x             - input signal
%     Fs            - sampling frequency (in Hz)
%     param_st      - parameter structure (with window length, window type, overlap, 
%                     frequency bands, total frequency band, and method)
%     SCALE_PSD     - scale PSD by factor of 2 (expect DC and Nyquist)? (0=no [default]
%                     or 1=yes)
%
% Outputs: 
%     pxx              - spectral estimate (e.g. power spectral density)
%     itotal_bandpass  - indices for total frequency band (defined in input args.)
%     f_scale          - frequency scaling factor
%     Nfreq            - total length of pxx (including -ve frequencies)
%     fp               - vector of frequencies (for plotting)
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%
%     param_st.method='robust-PSD';
%     param_st.L_window=2;
%     param_st.window_type='hamm';
%     param_st.overlap=50;
%     param_st.freq_bands=[0.5 4; 4 7; 7 13; 13 30];
%     param_st.total_freq_bands=[0.5 30];    
%
%     [pxx,itotal_bandpass,~,~,fp]=gen_spectrum(x,Fs,param_st);
%
%     figure(1); clf; hold all;
%     plot(fp,10*log10(pxx));
%     line([fp(itotal_bandpass(1)) fp(itotal_bandpass(1))],ylim,'color','k');        
%     line([fp(itotal_bandpass(end)) fp(itotal_bandpass(end))],ylim,'color','k');


% John M. O' Toole, University College Cork
% Started: 16-03-2017
%
% last update: Time-stamp: <2019-03-19 14:45:15 (otoolej)>
%-------------------------------------------------------------------------------
function [pxx,itotal_bandpass,f_scale,Nfreq,fp]=gen_spectrum(x,Fs,param_st,SCALE_PSD)
if(nargin<4 || isempty(SCALE_PSD)), SCALE_PSD=0; end

L_window=param_st.L_window;
window_type=param_st.window_type;
overlap=param_st.overlap;
spec_method=param_st.method;
    
    

% remove NaNs:
x(isnan(x))=[];

if(strcmp(lower(spec_method), 'bartlett-psd'))
    %---------------------------------------------------------------------
    % Bartlett PSD: same as Welch with 0% overlap and rectangular
    % window
    %---------------------------------------------------------------------
    window_type='rect';
    overlap=0;
    
    spec_method = 'psd';
end    


switch lower(spec_method)
  case 'psd'
    %---------------------------------------------------------------------
    % Welch PSD
    %---------------------------------------------------------------------
    [S_stft,Nfreq,f_scale,win_epoch]=gen_STFT(x,L_window,window_type,overlap,Fs);

    % average over time:
    pxx = nanmean(S_stft, 1)';
    N=length(pxx);

    % normalise (so similar to pwelch):
    E_win=sum(abs(win_epoch).^2)./Nfreq;
    pxx=(pxx./(Nfreq*E_win*Fs));

    
  case 'robust-psd'
    %---------------------------------------------------------------------
    % Welch PSD with median instead of mean
    %---------------------------------------------------------------------
    [S_stft,Nfreq,f_scale,win_epoch]=gen_STFT(x,L_window,window_type,overlap,Fs);

    % average over time:
    pxx = nanmedian(S_stft, 1)';    
    N=length(pxx);
    
    % normalise (so is similar to Welch's PSD):
    E_win=sum(abs(win_epoch).^2)./Nfreq;
    pxx=(pxx./(Nfreq*E_win*Fs));

    
    
  case 'periodogram'
    %---------------------------------------------------------------------
    % Periodogram
    %---------------------------------------------------------------------
    X=abs(fft(x)).^2;

    % +ve frequencies only:
    N=length(X); Nh=floor(N/2); Nfreq=N;
    X=X(1:Nh+1)';

    pxx=X./(Fs*N);


  otherwise
    fprintf('unknown spectral method ''%s''; check spelling\n',spec_method);
    pxx=NaN; itotal_bandpass=NaN; f_scale=NaN; fp=NaN;
end



% if need to scale (when calculating total power)
if(SCALE_PSD)
    pscale=ones(length(pxx),1)+1;     
    if(rem(Nfreq,2))
        pscale(1)=1;
    else
        pscale([1 end])=1;
    end
    pxx=pxx.*pscale;
end


N=length(pxx);
f_scale=(Nfreq/Fs);
% for plotting only:
fp=(0:(N-1))./f_scale;


if(isfield(param_st, 'total_freq_bands'))
    total_freq_bands=param_st.total_freq_bands;    

    % b) limit to frequency band of interest:
    itotal_bandpass=ceil(total_freq_bands(1)*f_scale):floor(total_freq_bands(2)*f_scale);
    itotal_bandpass=itotal_bandpass+1;
    itotal_bandpass(itotal_bandpass<1)=1;  itotal_bandpass(itotal_bandpass>N)=N;    
else
    itotal_bandpass=NaN;
end
