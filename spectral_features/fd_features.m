%-------------------------------------------------------------------------------
% fd_features: fractal dimension estimates
%
% Syntax: featx=fd_features(x,Fs,feat_name,params_st)
%
% Inputs: 
%     x,Fs,feat_name,params_st - 
%
% Outputs: 
%     featx - 
%
% Example:
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
% last update: Time-stamp: <2016-10-03 09:10:28 (otoolej)>
%-------------------------------------------------------------------------------
function featx=fd_features(x,Fs,feat_name,params_st)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='higuchi'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end

DBplot=0;

if(isempty(params_st))
    qEEGfs_parameters;
    if(strfind(feat_name,'fd'))
        params_st=feat_params_st.fd;
    else
        params_st=feat_params_st.(char(feat_name));
    end
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
        featx(n)=fd_higuchi(x,params_st.kmax);

      case 'katz'
        %---------------------------------------------------------------------
        % Katz estimate of fractal dimension [2]
        %---------------------------------------------------------------------
        x_diff=diff(x);
        L=sum( abs(x_diff) );
        d=max( abs(x(1)-x(2:end)) );
        p=length(x)-1;
        
        featx(n)=log10(p)/ (log10(d/L)+log10(p));

      otherwise
        error(['unknown feature: ' feat_name]);
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
        
        L(m)=sum( abs( x(m+ik.*k) - x(m+(ik-1).*k) ) )*(scale_factor/k);
        
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

    
