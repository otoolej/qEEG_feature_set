%-------------------------------------------------------------------------------
% channel_hemisphere_pairs: pair left/right channels
%
% Syntax: ipairs=channel_hemisphere_pairs(channel_labels)
%
% Inputs: 
%     channel_labels - 
%
% Outputs: 
%     ipairs - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 20-04-2016
%
% last update: Time-stamp: <2016-04-20 16:35:36 (otoolej)>
%-------------------------------------------------------------------------------
function ipairs=channel_hemisphere_pairs(channel_labels)

N=length(channel_labels);
[ileft,iright]=channel_hemispheres(channel_labels);

ipairs=NaN(2,length(ileft));

channel_labels=upper(channel_labels);


for n=1:length(ileft)
    ipairs(1,n)=ileft(n);    
    ch_left=upper(channel_labels{ileft(n)});
    
    % change numbers from odd to even
    ch_left_match= ...
        strrep(strrep(strrep(strrep(ch_left,'1','2'),'3','4'),'5','6'),'7','8');
    
    imatch=find(strcmp(channel_labels(iright),ch_left_match));

    % and check for reversed order:
    sep=strfind(ch_left_match,'-');% sep=sep{1};
    ch1=ch_left_match(1:(sep-1)); ch2=ch_left_match((sep+1):end);
    ch_left_match_rv=[ch2 '-' ch1];
    
    imatch_rv=find(strcmp(channel_labels(iright),ch_left_match_rv));
    
    if(~isempty(imatch))
        ipairs(2,n)=iright(imatch);        
    elseif(~isempty(imatch_rv))
        ipairs(2,n)=iright(imatch_rv);        
    else
        warning(sprintf('no matching pair for channel: %s',ch_left));
    end


% $$$     dispVars(channel_labels(ipairs(1,n)),channel_labels(ipairs(2,n)));
end


