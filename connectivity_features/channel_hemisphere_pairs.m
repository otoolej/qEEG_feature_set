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
% last update: Time-stamp: <2016-04-27 10:23:32 (otoolej)>
%-------------------------------------------------------------------------------
function ipairs=channel_hemisphere_pairs(channel_labels)

DBverbose=0;

N=length(channel_labels);
[ileft,iright]=channel_hemispheres(channel_labels);

N_left=length(ileft);

ipairs=NaN(2,N_left);

channel_labels=upper(channel_labels);


for n=1:N_left
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
        ipairs(1:2,n)=NaN;
        if(DBverbose)
            fprintf('no matching pair for channel: %s\n',ch_left);        
% $$$         warning(sprintf('no matching pair for channel: %s',ch_left));
        end
    end


end

irem=[];
for n=1:N_left
    if(any(isnan(ipairs(:,n))))
        irem=[irem n];
    end
end
if(~isempty(irem))
    ipairs(:,irem)=[];    
    
    if(DBverbose),  print_table(channel_labels(ipairs)); end
end





