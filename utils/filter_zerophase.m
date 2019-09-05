%-------------------------------------------------------------------------------
% filter_zerophase: Simple implementation of zero-phase FIR filter using 'filtfilt'
%
% Syntax: x_filt=filter_zerophase(x,Fs,LP_fc,HP_fc,L_filt)
%
% Inputs:
%     x      - input signal 
%     Fs     - sample frequency (Hz)
%     LP_fc  - lowpass cut off (Hz)
%     HP_fc  - highpass cut off (Hz)
%     L_filt - length of filter (in samples)
%
% Outputs: 
%     y - filtered signal
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%     LP_fc=20; HP_fc=1;
%
%     y=filter_zerophase(x,Fs,LP_fc,HP_fc,501);
%
%     figure(1); clf; hold all;
%     ttime=(0:(length(x)-1))./Fs;
%     plot(ttime,x); plot(ttime,y);
%

% John M. O' Toole, University College Cork
% Started: 08-05-2013
%
% last update: Time-stamp: <2019-09-05 14:03:09 (otoolej)>
%-------------------------------------------------------------------------------
function x_filt=filter_zerophase(x,Fs,LP_fc,HP_fc,L_filt,win_type,DBplot)
if(nargin<2 || isempty(Fs)), Fs=1; end
if(nargin<3 || isempty(LP_fc)), LP_fc=[]; end
if(nargin<4 || isempty(HP_fc)), HP_fc=[]; end
if(nargin<5 || isempty(L_filt)), L_filt=[]; end
if(nargin<6 || isempty(win_type)), win_type=[]; end
if(nargin<7 || isempty(DBplot)), DBplot=0; end


if(LP_fc==0), LP_fc=[]; end
if(HP_fc==0), HP_fc=[]; end


% replace NaNs with zeros:
inans=find(isnan(x));
if(~isempty(inans))
    x(inans)=0;
end


N=length(x);
if(isempty(L_filt) || L_filt>(N/4))
    L_filt=make_odd(round(N/4-2)); 
end

wwin=[];
if(~isempty(win_type))
    wwin=gen_window(win_type,L_filt);
    
    fir_ = @(L, f_cut, ftype, wwin) fir1(L, f_cut, ftype, wwin);
else
    % update for R2019a: fir1.m no-long ignores empty window:
    fir_ = @(L, f_cut, ftype, wwin) fir1(L, f_cut, ftype);
end

%---------------------------------------------------------------------
% either bandpass, low-pass, or high-pass
%---------------------------------------------------------------------
if(~isempty(LP_fc) && ~isempty(HP_fc))
    
    if(LP_fc>HP_fc)
        % band-pass:
        filt_passband=[HP_fc/(Fs/2) LP_fc/(Fs/2)];
        b = fir_(L_filt-1, filt_passband, 'bandpass', wwin);  
    else
        % or band-stop:
        filt_passband=[LP_fc/(Fs/2) HP_fc/(Fs/2)];
        b = fir_(L_filt-1, filt_passband, 'stop', wwin);  
    end
        
        
elseif(~isempty(LP_fc))
    filt_passband=LP_fc/(Fs/2);
    b = fir_(L_filt-1, filt_passband, 'low', wwin);

elseif(~isempty(HP_fc))
    L_filt=L_filt+1;
    filt_passband=HP_fc/(Fs/2);
    b = fir_(L_filt-1, filt_passband, 'high', wwin);
    
else
    error('need to specify cut-off frequency.');
end

  
% Do the filtering
% $$$ x_filt=filtfilt(b,1,[x(:); zeros(L_filt,1)]);
x_filt=filtfilt(b,1,x(:));

x_filt=x_filt(1:N);

if(DBplot)
    fvtool(b,1,'Fs',Fs);
    
    figure(10); clf; 
    plot(x,'b'); hold on;
    plot(x_filt,'r');
    
    figure(11); clf; hold all;
    pfftva(x,'Fs',Fs,'db',1);    
    pfftva(x_filt,'Fs',Fs,'db',1);
end


% swap back in the NaNs:
if(~isempty(inans))
    x_filt(inans)=NaN;
end



function w=gen_window(win_type,L)
%---------------------------------------------------------------------
% window
%---------------------------------------------------------------------
win_type=lower(win_type);

switch win_type
  case {'hamm','hamming'}
    w=hamming(L);
  case {'hann','hanning'}
    w=hanning(L);
  case {'rect','rectangular'}
    w=ones(1,L);
  case {'black','blackman'}
    w=blackman(L);
  case {'blackmanharris'}
    w=blackmanharris(L);
  case {'kaiser'}
    w=kaiser(L,5);
  otherwise
    wfn=str2func(win_type);
    w=wfn(L);
end

    
