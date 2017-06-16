%-------------------------------------------------------------------------------
% filter_butterworth_withnans: Filter using butterworth filter with zero-phase filter
%
% Syntax: y=filter_butterworth_withnans(x,Fs,F3db_lowpass,F3db_highpass)
%
% Inputs: 
%     x             - input signal 
%     Fs            - sample frequency (Hz)
%     F3db_lowpass  - 3 dB lowpass cut off (Hz)
%     F3db_highpass - 3 dB highpass cut off (Hz)
%     order         - filter order
%     FILTER_REPLACE_ARTEFACTS - what to do with NaNs?
%         either replace with 0s ('zeros'), linear interpolation 
%         ('linear_interp', default), or cubic spline interpolation ('cubic_interp')
%
% Outputs: 
%     y - filtered signal
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%     F3db_lowpass=30; F3db_highpass=0.5;
%
%     y=filter_butterworth_withnans(x,Fs,F3db_lowpass,F3db_highpass,5);
%
%     figure(1); clf; hold all;
%     ttime=(0:(length(x)-1))./Fs;
%     plot(ttime,x); plot(ttime,y);


% John M. O' Toole, University College Cork
% Started: 31-10-2013
%
% last update: Time-stamp: <2017-06-15 17:47:21 (otoolej)>
%-------------------------------------------------------------------------------
function [y,inans]=filter_butterworth_withnans(x,Fs,F3db_lowpass,F3db_highpass, ...
                                       order,FILTER_REPLACE_ARTEFACTS,DB)
if(nargin<2 || isempty(Fs)), Fs=1; end
if(nargin<3 || isempty(F3db_lowpass) || F3db_lowpass==0), F3db_lowpass=[]; end
if(nargin<4 || isempty(F3db_highpass) || F3db_highpass==0), F3db_highpass=[]; end
if(nargin<5 || isempty(order)), order=3; end
if(nargin<6 || isempty(FILTER_REPLACE_ARTEFACTS)), FILTER_REPLACE_ARTEFACTS='linear_interp'; end
if(nargin<7 || isempty(DB)), DB=0; end

inans=[];

% check values:
if(F3db_lowpass<=0 | F3db_lowpass>=(Fs/2) | ...
   F3db_highpass<=0 | F3db_highpass>=(Fs/2))
    fprintf('invalid filter value; ignoring\n');
    fprintf('  (high-pass/low-pass cut off: %g/%g)\n',F3db_highpass,F3db_lowpass);
    y=x;
    return;
end


if(isempty(F3db_highpass))
    [b,a]=butter(order,F3db_lowpass/(Fs/2),'low');    
    
elseif(isempty(F3db_lowpass))
    [b,a]=butter(order,F3db_highpass/(Fs/2),'high');        
    
else
    if(length(order)>1)
        order_low=order(1);                
        order_high=order(2);
    else
        order_low=order;        
        order_high=order;        
    end
    
    [y,inans_low]=filter_butterworth_withnans(x,Fs,F3db_lowpass,[],order_low, ...
                                              FILTER_REPLACE_ARTEFACTS, DB);
    [y,inans_high]=filter_butterworth_withnans(y,Fs,[],F3db_highpass,order_high, ...
                                               FILTER_REPLACE_ARTEFACTS, DB);
    
    inans=unique([inans_low inans_low]);
    return;    
end


% remove NaNs and replace with ?
inans=find(isnan(x));
if(~isempty(inans))
    switch FILTER_REPLACE_ARTEFACTS
      case 'zeros'
        x(inans)=0;
        
      case 'linear_interp'
        x=replace_start_ends_NaNs_with_zeros(x);
        x=naninterp(x,'linear');
        
      case {'cubic_interp','nans'}
        x=replace_start_ends_NaNs_with_zeros(x);
        x=naninterp(x,'pchip');
    end
end



y=filtfilt(b,a,x);


% special case: if NaNs
if(strcmp(FILTER_REPLACE_ARTEFACTS,'nans'))
    if(~isempty(inans))
        y(inans)=NaN;
    end
end

if(DB)
    fvtool(b,a,'Fs',Fs); %,'frequencyscale','log');
    figure(9); clf; hold all;
    plot(1:length(x),x,1:length(x),y);
    legend({'original','filtered'});
end





function x=replace_start_ends_NaNs_with_zeros(x)
%---------------------------------------------------------------------
% replace leading or trailing NaNs with zeros (needed for naninterp.m)
%---------------------------------------------------------------------
N=length(x);

istart=find(~isnan(x),1,'first');
iend=find(~isnan(x),1,'last');

if(~isempty(istart) & istart>1)
    x(1:istart)=0;
end
if(~isempty(iend) & iend<N)
    x((iend+1):N)=0;
end


function [X,inan]=naninterp(X,method)
%---------------------------------------------------------------------
% fill 'gaps' in data (marked by NaN) by interpolating
%---------------------------------------------------------------------
if(nargin<2 || isempty(method)), method='linear'; end

inan=find(isnan(X));
if(isempty(inan))
  return;
elseif(length(inan)==1)
  if(inan>1)
    X(inan)=X(inan-1);
  else
    X(inan)=X(inan+1);
  end
else
    try
        X(inan)=interp1(find(~isnan(X)), X(~isnan(X)), inan, method);
    catch
        error('linear interpolation with NaNs');
    end
end

