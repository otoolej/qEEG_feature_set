%--------------------------------------------------------------------------------
% General function to calculate a window function (mostly a wrapper for the existing
% window functions).
% 
% USE: win = get_window( win_length, win_type, [win_param], [DFT_WINDOW], [Npad] )
%
% INPUT:
%       win_length = length of window
%       win_type   = type of window: { 'delta' | 'rect' | 'bart' | 'hamm'
%                    | 'hann' | 'tukey' | 'gauss' | 'cosh' }
%       win_param  = (optional) window parameter.
%       DFT_WINDOW = (optional) parameter { 0 | 1 }. 
%                     If 1 returns DFT of window. 
%       Npad       = (optional) zero-pad window to length Npad.
%
% OUTPUT:
%       win = window 
%
% EXAMPLE
%      N=64; 
%      win=get_window(N,'tukey',0.4,0,2*N); 
%      plot(win);


% Copyright (C) 2007,2008 John M. O' Toole, The University of Queensland
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------------
function win = get_window( win_length, win_type, win_param, DFT_WINDOW, Npad )
if( nargin<3 ) win_param=[]; end
if( nargin<4 ) DFT_WINDOW=0; end
if( nargin<5 ) Npad=0; end


win=get_win(win_length,win_type,win_param,DFT_WINDOW);
win=shift_win(win);

if(Npad>0)
  win=pad_win(win,Npad);
end



function win=get_win(win_length, win_type, win_param, DFT_WINDOW )
%--------------------------------------------------------------------------------
% Get the window. Negative indices are first.
%--------------------------------------------------------------------------------

switch win_type
 case {'delt', 'delta'}
  win=zeros(win_length,1);
  wh = floor(win_length/2);
  win(wh+1)=1;
 case {'rect' , 'rectangular'}
  win(1:win_length) = 1;
 case {'bart', 'bartlett'}
  win = bartlett( win_length );
 case {'hamm', 'hamming'}
  win = hamming( win_length );
 case {'hann', 'hanning'}
  win = hanning( win_length );
 case {'tuke', 'tukey' }
  % NOTE: seems to be problem with Octave's (v2.9.12) tukeywin.m for N odd.
  if(isempty(win_param))
    win = tukeywin( win_length );
  else
    win = tukeywin( win_length, win_param );
  end
 case {'gaus','gauss'}
  if(isempty(win_param))
    win = gausswin( win_length );
  else
    win = gausswin( win_length, win_param );
  end
 case 'cosh'
  win_hlf = fix( win_length / 2);

  if(isempty(win_param))
    win_param = 0.01;
  end
  for m = -win_hlf:win_hlf
    win(mod(m,win_length)+1) = cosh( m ).^( -2 * win_param );
  end
  win = fftshift(win);
  
  
 otherwise
  error(['Unknown window type ' win_type]);
end
  

%---------------------------------------------------------------------
% If want the DFT of win
%---------------------------------------------------------------------
if(DFT_WINDOW)
  win=circshift(win(:),ceil(win_length/2));
  win=fft(win);
  win=circshift(win(:),floor(win_length/2));
end


function w=shift_win(w)  
%--------------------------------------------------------------------------------
% Shift the window so that positive indices are first.
%--------------------------------------------------------------------------------
N=length(w);
w=circshift(w(:),ceil(N/2));



function w_pad=pad_win(w,Npad)
%--------------------------------------------------------------------------------
%
% Pad window to Npad. 
%
% Presume that positive window indices are first. 
%
% When N is even use method described in [1]
%
%   References:
%     [1] S. Lawrence Marple, Jr., Computing the discrete-time analytic
%     signal via FFT, IEEE Transactions on Signal Processing, Vol. 47,
%     No. 9, September 1999, pp.2600--2603.
%
%--------------------------------------------------------------------------------
w=w(:);
w_pad=zeros(Npad,1);
N=length(w);
Nh=floor(N/2);
if(Npad<N) error('Npad is less than N'); end


% Trival case:
if(N==Npad) 
  w_pad=w; 
  return; 
end


 % For N odd:
if( rem(N,2)==1 )
  n=0:Nh;
  w_pad(n+1)=w(n+1); 
  n=1:Nh;
  w_pad(Npad-n+1)=w(N-n+1);
 
  % For N even:
  % split the Nyquist frequency in two and distribute over positive
  % and negative indices.
else
  n=0:(Nh-1);
  w_pad(n+1)=w(n+1); 
  w_pad(Nh+1)=w(Nh+1)/2; 

  n=1:Nh-1;
  w_pad(Npad-n+1)=w(N-n+1);
  w_pad(Npad-Nh+1)=w(Nh+1)/2; 
end
