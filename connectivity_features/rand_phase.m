%-------------------------------------------------------------------------------
% rand_phase: randomize phase (of FFT) for real-valued signal
%
% Syntax: y=rand_phase(x,N_iter)
%
% Inputs: 
%     x      - input signal, length-N
%     N_iter - number of random signals required
%
% Outputs: 
%     y - output signal, same amplitude spectrum but with random phase 
%         (size N_iter x N)
%
% Example:
%     N=1000; N_iter=500;
%     x=randn(1,N);
%     y=rand_phase(x,N_iter);
%
%     figure(1); clf; hold all;
%     subplot(2,1,1);
%     plot(x);
%     subplot(2,1,2);
%     plot(y);
%
% 

% John M. O' Toole, University College Cork
% Started: 24-08-2016
%
% last update: Time-stamp: <2019-09-05 13:59:58 (otoolej)>
%-------------------------------------------------------------------------------
function y=rand_phase(x,N_iter)
if(nargin<2 || isempty(N_iter)), N_iter=1; end


if(size(x,2)==1), x=x'; end

[L,N]=size(x);

N=length(x);
Nh=floor(N/2);


X=fft(x);

%---------------------------------------------------------------------
% 1. generate random phase
%---------------------------------------------------------------------
if(rem(N,2))
    rphase=-pi+2*pi*rand(N_iter,Nh); 
    
    rphase=[zeros(N_iter,1) rphase(:,1:Nh) -rphase(:,end:-1:1)];
else
    rphase=-pi+2*pi*rand(N_iter,Nh-1); 
    
    rphase=[zeros(N_iter,1) rphase zeros(N_iter,1) -rphase(:,end:-1:1)];
end
        

%---------------------------------------------------------------------
% apply to spectrum and FFT back to time domain:
%---------------------------------------------------------------------
if(N_iter>1)
    Y=bsxfun(@times,exp(j.*rphase),abs(X));
else
    Y=abs(X).*exp(j.*rphase);
end

y=real( ifft(Y.') )';


DBplot=0;
if(DBplot)
    figure(18); clf; hold all;
    plot(x);
    plot(y');

    if(N_iter>1), yplot=y(end,:); end
    
    figure(19); clf; hold all;
    subplot(2,2,[1 2]); hold all;
    plot(x); plot(yplot);

    subplot(2,2,3); hold all;
    plot(abs(X)); plot(abs(fft(yplot)));
    subplot(2,2,4); hold all;
    plot(angle(X)); plot(angle(fft(yplot)));
    
    figure(9); clf; hold all;
    hists_two(x,yplot);
    
    dispEE(abs(fft(yplot)),abs(X));
end

