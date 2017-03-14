%---------------------------------------------------------------------
% run this to add paths
%---------------------------------------------------------------------

% what directory is this file in?
cur_full_path=mfilename('fullpath');
cur_fname=mfilename();
i=findstr(cur_full_path,cur_fname);
cur_dir=cur_full_path(1:(i-1));

% add the paths:
addpath( genpath(cur_dir) );


clear cur_full_path cur_fname cur_dir;

