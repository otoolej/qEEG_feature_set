%-------------------------------------------------------------------------------
% len_cont_zeros: find length of continuous segments of zeros from binary mask x. Can
% contain NaNs.
%
%
% Syntax: [lens,istart,iend]=len_cont_zeros(x)
%
% Inputs: 
%     x     - binary [0,1] vector
%     const - which to look for, either 0 (default) or 1
%
% Outputs: 
%     [lens,istart,iend] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 25-02-2014
%
% last update: Time-stamp: <2014-02-28 18:06:55 (otoolej)>
%-------------------------------------------------------------------------------
function [lens,istart,iend]=len_cont_zeros(x,const)
if(nargin<2 || isempty(const)), const=0; end

DBplot=0;

x=x(:).';

if( ~all(ismember(sort(unique(x(~isnan(x)))),[0 1])) || ...
    ~ismember(const,[0 1]) )
    warning('must be binary signal');
    return;
end
if(const==1)
    y=invert_bin_array_with_NaNs(x); 
else
    y=x;
end

% find run of zeros:
iedge=diff([0 y==0 0]);
istart=find(iedge==1);
iend=find(iedge==-1)-1;
lens=[iend-istart];


if(DBplot)
    figure(100); clf; hold all;
    plot(x);  
    plot(istart,x(istart),'x');
    plot(iend,x(iend),'o');    
    ylim([-.1 1.2]);
end



function y=invert_bin_array_with_NaNs(x)
y=x;
y(~isnan(x))=~x(~isnan(x));
