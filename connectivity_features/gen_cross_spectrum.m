%-------------------------------------------------------------------------------
% gen_cross_spectrum: cross-spectrums
%
% Syntax: [pxy,Nfreq,f_scale,fp]=gen_cross_spectrum(x,y,Fs,param_st)
%
% Inputs: 
%     x,y,Fs,param_st - 
%
% Outputs: 
%     [pxy,Nfreq,f_scale,fp] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 17-03-2017
%
% last update: Time-stamp: <2017-03-28 17:21:02 (otoolej)>
%-------------------------------------------------------------------------------
function [pxy,Nfreq,f_scale,fp]=gen_cross_spectrum(x,y,Fs,param_st)

L_window=param_st.L_window;
window_type=param_st.window_type;
overlap=param_st.overlap;
freq_bands=param_st.freq_bands;
spec_method=param_st.method;
    
    

% remove NaNs:
x(isnan(x))=[];


switch lower(spec_method)
  case 'psd'
    %---------------------------------------------------------------------
    % Welch cross-PSD 
    %---------------------------------------------------------------------
    [S_x,Nfreq,f_scale,win_epoch]=gen_STFT(x,L_window,window_type,overlap,Fs,1);
    [S_y,Nfreq,f_scale,win_epoch]=gen_STFT(y,L_window,window_type,overlap,Fs,1);

    S_xy=S_x.*conj(S_y);

    pxy=nanmean(S_xy)';
    N=length(pxy);

    % normalise (so is similar to Welch's PSD):
    E_win=sum(abs(win_epoch).^2)./Nfreq;
    pxy=(pxy./(Nfreq*E_win*Fs));
    
    
  case 'robust-psd'
    %---------------------------------------------------------------------
    % similar to Welch but use median instead of mean value
    %---------------------------------------------------------------------
    [S_x,Nfreq,f_scale,win_epoch]=gen_STFT(x,L_window,window_type,overlap,Fs,1);
    [S_y,Nfreq,f_scale,win_epoch]=gen_STFT(y,L_window,window_type,overlap,Fs,1);

    S_xy=S_x.*conj(S_y);

    pxy=nanmedian(S_xy)';
    N=length(pxy);
    
    % normalise (so is similar to Welch's PSD):
    E_win=sum(abs(win_epoch).^2)./Nfreq;
    pxy=(pxy./(Nfreq*E_win*Fs));

    
  case 'periodogram'
    %---------------------------------------------------------------------
    % Periodogram
    %---------------------------------------------------------------------
    X=fft(x);
    Y=fft(y);    
    pxy=X.*conj(Y);

    % +ve frequencies only:
    N=length(pxy); Nh=floor(N/2); Nfreq=N;
    pxy=pxy(1:Nh+1)';

    pxy=pxy./(Fs*N);
    
    
  otherwise
    fprintf('unknown spectral method ''%s''; check spelling\n',spec_method);
    pxy=NaN; f_scale=NaN; fp=NaN;
end


N=length(pxy);
f_scale=(Nfreq/Fs);
% for plotting only:
fp=(0:(N-1))./f_scale;

