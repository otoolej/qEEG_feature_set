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
% last update: Time-stamp: <2019-03-19 15:03:23 (otoolej)>
%-------------------------------------------------------------------------------
function [pxy,Nfreq,f_scale,fp]=gen_cross_spectrum(x,y,Fs,param_st)

L_window=param_st.L_window;
window_type=param_st.window_type;
overlap=param_st.overlap;
freq_bands=param_st.freq_bands;
spec_method=param_st.method;
    
    
% remove NaNs:
x(isnan(x))=[];

if(strcmp(lower(spec_method), 'bartlett-psd'))
    %---------------------------------------------------------------------
    % Bartlett PSD: same as Welch with 0% overlap and rectangular
    % window
    %---------------------------------------------------------------------
    window_type = 'rect';
    overlap = 0;
    
    spec_method = 'psd';
end    




%---------------------------------------------------------------------
% Welch cross-PSD 
%---------------------------------------------------------------------
[S_x,Nfreq,f_scale,win_epoch]=gen_STFT(x,L_window,window_type,overlap,Fs,1);
[S_y,Nfreq,f_scale,win_epoch]=gen_STFT(y,L_window,window_type,overlap,Fs,1);

S_xy=S_x.*conj(S_y);


switch lower(spec_method)
  case 'psd'
    %---------------------------------------------------------------------
    % mean for Welch PSD
    %---------------------------------------------------------------------
    pxy = nanmean(S_xy, 1).';
    
  case 'robust-psd'
    %---------------------------------------------------------------------
    % median for robust PSD
    %---------------------------------------------------------------------
    pxy = nanmedian(S_xy, 1).';
    
  otherwise
    fprintf('unknown cross-spectral method ''%s''; check spelling\n',spec_method);
    pxy=NaN; f_scale=NaN; fp=NaN;
end


% normalise (so is similar to Welch's PSD):
E_win=sum(abs(win_epoch).^2)./Nfreq;
pxy=(pxy./(Nfreq*E_win*Fs));


N=length(pxy);
f_scale=(Nfreq/Fs);
% for plotting only:
fp=(0:(N-1))./f_scale;

