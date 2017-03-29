%-------------------------------------------------------------------------------
% gen_STFT: Short-time Fourier transform (or spectrogram)
%
% Syntax: [S_stft,Nfreq,f_scale,win_epoch]=gen_STFT(x,L_window,window_type,overlap,Fs)
%
% Inputs: 
%     x            - input signal
%     L_window     - window length
%     window_type  - window type
%     overlap      - percentage overlap
%     Fs           - sampling frequency (Hz)
%     STFT_OR_SPEC - return short-time Fourier transform (STFT) or spectrogram
%                    (0=spectrogram [default] and 1=STFT)
%
% Outputs: 
%     S_stft     - spectrogram 
%     Nfreq      - length of FFT
%     f_scale    - frequency scaling factor
%     win_epoch  - window
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%
%     L_window=2;
%     window_type='hamm';
%     overlap=80;
%
%     S_stft=gen_STFT(x,L_window,window_type,overlap,Fs);
%
%     figure(1); clf; hold all;
%     imagesc(S_stft); axis('tight');



% John M. O' Toole, University College Cork
% Started: 16-03-2017
%
% last update: Time-stamp: <2017-03-28 18:22:47 (otoolej)>
%-------------------------------------------------------------------------------
function [S_stft,Nfreq,f_scale,win_epoch]=gen_STFT(x,L_window,window_type,overlap,Fs, ...
                                                  STFT_OR_SPEC)
if(nargin<6 || isempty(STFT_OR_SPEC)), STFT_OR_SPEC=0; end


[L_hop,L_epoch,win_epoch]=gen_epoch_window(overlap,L_window,window_type,Fs,1);


N=length(x);
N_epochs=floor( (N-(L_epoch-L_hop))/L_hop );
if(N_epochs<1) N_epochs=1; end
nw=0:L_epoch-1;
Nfreq=L_epoch;


%---------------------------------------------------------------------
% generate short-time FT on all data:
%---------------------------------------------------------------------
K_stft=zeros(N_epochs,L_epoch);
for k=1:N_epochs
    nf=mod(nw+(k-1)*L_hop,N);
    
    K_stft(k,:)=x(nf+1).*win_epoch';
end

f_scale=(Nfreq/Fs);
if(STFT_OR_SPEC)
    S_stft=fft(K_stft.',Nfreq);    
else
    S_stft=abs(fft(K_stft.',Nfreq)).^2;
end

S_stft=S_stft';

% return only positive frequencies:
S_stft=S_stft(:,1:floor(Nfreq/2)+1);










