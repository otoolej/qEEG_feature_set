%-------------------------------------------------------------------------------
% filter_butterworth_withnans: Filter using butterworth filter with zero-phase filter
%
% Syntax: y=filter_butterworth_withnans(x,Fs,F3db_lowpass,F3db_highpass)
%
% Inputs: 
%     x,Fs,F3db_lowpass,F3db_highpass - 
%
% Outputs: 
%     y - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 31-10-2013
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
    [y,inans_low]=filter_butterworth_withnans(x,Fs,F3db_lowpass,[],order, ...
                                              FILTER_REPLACE_ARTEFACTS, DB);
    [y,inans_high]=filter_butterworth_withnans(y,Fs,[],F3db_highpass,order, ...
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
% replace leading or trailing NaNs with zeros (needed to naninterp.m)
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
