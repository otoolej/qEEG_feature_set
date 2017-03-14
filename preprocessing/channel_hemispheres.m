%-------------------------------------------------------------------------------
% channel_hemispheres: return indices for left and right channels
%
% Syntax: [ileft,iright]=channel_hemispheres(channels_all)
%
% Inputs: 
%     channel_labels - cell of bipolar channel names
%                      e.g. {'C3-O1','C4-O2', 'F3-C3', 'F4-C4'}
%
% Outputs: 
%     ileft  - indices of the channels on left hemispheres
%     iright - indices of the channels on right hemispheres
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     channel_labels=data_st.ch_labels;
%     
%     [ileft,iright]=channel_hemispheres(channel_labels);
%
%     fprintf('left hemisphere channels: %s\n',strjoin(channel_labels(ileft),', '));
%     fprintf('right hemisphere channels: %s\n',strjoin(channel_labels(iright),', '));
%

% John M. O' Toole, University College Cork
% Started: 02-03-2015
%
% last update: Time-stamp: <2017-03-14 14:52:51 (otoolej)>
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

