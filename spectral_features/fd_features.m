%-------------------------------------------------------------------------------
% fd_features: fractal dimension estimates
%
% Syntax: featx=fd_features(x,Fs,params_st)
%
% Inputs: 
%     x          - epoch of EEG data (size 1 x N)
%     Fs         - sampling frequency (in Hz)
%     params_st  - parameters (as structure); see neural_parameters.m for examples
%
% Outputs: 
%     featx  - fractal dimension for each frequency band 
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%
%     featx=fd_features(x,Fs);
%
%     
%
% [1] T Higuchi, “Approach to an irregular time series on the basis of the fractal
% theory,” Phys. D Nonlinear Phenom., vol. 31, pp. 277–283, 1988.
%
% [2] MJ Katz, Fractals and the analysis of waveforms. Computers in Biology and Medicine,
% vol. 18, no. 3, pp. 145–156. 1988


% John M. O' Toole, University College Cork
% Started: 03-10-2016
%
% last update: Time-stamp: <2018-09-18 10:56:25 (otoolej)>
%-------------------------------------------------------------------------------
function featx=fd_features(x,Fs,params_st)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(params_st)), params_st=[]; end

DBplot=0;

if(isempty(params_st))
    neural_parameters;
    params_st=feat_params_st.FD;
end

freq_bands=params_st.freq_bands;
N_freq_bands=size(freq_bands,1);
if(isempty(freq_bands))
    N_freq_bands=1;
end

x_orig=x;

for n=1:N_freq_bands
    if(~isempty(freq_bands))
        x=filter_butterworth_withnans(x_orig,Fs,freq_bands(n,2),freq_bands(n,1),5, ...
                                      params_st.FILTER_REPLACE_ARTEFACTS);
    end

    switch params_st.method
      case 'higuchi'
        %---------------------------------------------------------------------
        % Higuchi estimate of fractal dimension [1]
        %---------------------------------------------------------------------
        featx(n)=fd_higuchi(x,params_st.qmax);

      case 'katz'
        %---------------------------------------------------------------------
        % Katz estimate of fractal dimension [2]
        %---------------------------------------------------------------------
        featx(n)=fd_katz(x);

      otherwise
        fprintf('unknown feature ''%s''; check spelling\n',feat_name);
        featx=NaN;
    end

end


function [FD,r2,k_all,L_avg]=fd_higuchi(x,kmax)
%---------------------------------------------------------------------
% Higuchi estimate in [1]
%---------------------------------------------------------------------
N=length(x);

DBplot=0;

if(nargin<2 || isempty(kmax)), kmax=floor(N/10); end

FD=[]; L=[];


% what values of k to compute?
ik=1; k_all=[]; knew=0;
while( knew<kmax )
    if(ik<=4)
        knew=ik;
    else
        knew=floor(2^((ik+5)/4));
    end
    if(knew<=kmax)
        k_all=[k_all knew];
    end
    ik=ik+1;
end


%---------------------------------------------------------------------
% curve length for each vector:
%---------------------------------------------------------------------
inext=1; L_avg=zeros(1,length(k_all));
for k=k_all
    
    L=zeros(1,k);
    for m=1:k
        ik=1:floor( (N-m)/k );
        scale_factor=(N-1)/(floor( (N-m)/k )*k);
        
        L(m)=nansum( abs( x(m+ik.*k) - x(m+(ik-1).*k) ) )*(scale_factor/k);
    end

    L_avg(inext)=mean(L);
    inext=inext+1;    
end

x1=log2(k_all); y1=log2(L_avg);
c=polyfit(x1,y1,1);
FD=-c(1);



if(nargout>1)
    y_fit=c(1)*x1 + c(2);
    y_residuals=y1-y_fit;

    r2=1-(sum( y_residuals.^2 ))./( (N-1).*var(y1) );
end


if(DBplot)
    figure(42); clf; hold all;
    plot(log2(k_all),log2(L_avg),'o');
    y_fit=c(1)*log2(k_all) + c(2);
    plot(log2(k_all),y_fit,'-');
% $$$     ylim([-5 15]);
end

    


function D=fd_katz(x)
%---------------------------------------------------------------------
% Katz estimate in [2]
%---------------------------------------------------------------------
N=length(x);
p=N-1;

% 1. line-length
for n=1:N-1
    L(n)=sqrt( 1 + (x(n)-x(n+1)).^2 );
end
L=sum(L);

% 2. maximum distance:
d=zeros(1,p);
for n=1:N-1
    d(n)=sqrt( n.^2 + (x(1)-x(n+1)).^2 );
end
d=max(d);



D=log(p)/(log(d/L)+log(p));
