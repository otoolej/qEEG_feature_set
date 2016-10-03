%-------------------------------------------------------------------------------
% get_epoch_window: calculate overlap size (in samples) and window length for
% overlap-and-add type analysis
%
% Syntax: [L_hop,L_epoch,win_epoch]=get_epoch_window(L_overlap,L_epoch,win_type,Fs)
%
% Inputs: 
%     L_overlap,L_epoch,win_type,Fs - 
%
% Outputs: 
%     [L_hop,L_epoch,win_epoch] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 28-05-2013
%-------------------------------------------------------------------------------
function [L_hop,L_epoch,win_epoch]=get_epoch_window(L_overlap,L_epoch,win_type,Fs)
L_hop=(100-L_overlap)/100;


L_epoch=floor( L_epoch*Fs );

% check for window type to force constant-overlap add constraint
% i.e. \sum_m w(n-mR)=1 for all n, where R=overlap size
win_type=lower(win_type(1:4));


switch win_type
  case 'hamm'
    L_hop=(L_epoch-1)*L_hop;
  case 'hann'
    L_hop=(L_epoch+1)*L_hop;
  otherwise
    L_hop=L_epoch*L_hop;
end

L_hop=ceil(L_hop);

win_epoch=get_window(L_epoch,win_type);
win_epoch=circshift(win_epoch,floor(length(win_epoch)/2));



% special case to force constant overlap-add constraint:
% (see SPECTRAL AUDIO SIGNAL PROCESSING, JULIUS O. SMITH III)
if( strcmpi(win_type,'hamm')==1 && rem(L_epoch,2)==1 )
    win_epoch(1)=win_epoch(1)/2;
    win_epoch(end)=win_epoch(end)/2;
end    

if( strcmpi(win_type,'rect')==1 ) win_epoch=win_epoch.'; end
