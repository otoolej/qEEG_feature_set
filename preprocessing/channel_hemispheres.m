%-------------------------------------------------------------------------------
% channel_hemispheres: return indices for left and right channels
%
% Syntax: [ileft,iright]=channel_hemispheres(channels_all)
%
% Inputs: 
%     channels_all - 
%
% Outputs: 
%     [ileft,iright] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 02-03-2015
%
% last update: Time-stamp: <2017-03-13 11:10:30 (otoolej)>
%-------------------------------------------------------------------------------
function [ileft,iright]=channel_hemispheres(channels_all)


N_channels=length(channels_all);

ileft=[];  iright=[];
for n=1:N_channels
    cname=channels_all{n};
    M=length(cname);
    il=[]; ir=[];
    for p=1:M
        num=str2num(cname(p));
        if(~isempty(num))
            if(rem(num,2))
                il=[il n];
            else
                ir=[ir n;];
            end
        end
    end
    
    if(~isempty(il) && ~isempty(ir))
        error(['both odd and even in channel: ' cname]);
    end
    if(~isempty(il))
        ileft=[ileft n];
    elseif(~isempty(ir))        
        iright=[iright n];
    else
        warning(['left or right channel: ' cname '?']);
    end
end

