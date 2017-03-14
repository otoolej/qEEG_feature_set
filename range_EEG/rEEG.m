%-------------------------------------------------------------------------------
% rEEG: range EEG as defined in [1]
%
% Syntax: featx=rEEG(x,Fs,feat_name,params_st)
%
% Inputs: 
%     x          - epoch of EEG data (size 1 x N)
%     Fs         - sampling frequency (in Hz)
%     feat_name  - feature type, defaults to 'rEEG_mean';
%                  see full list of 'rEEG_' features in all_features_list.m
%     params_st  - parameters (as structure); 
%                  see neural_parameters.m for examples
%
% Outputs: 
%     featx  - feature at each frequency band 
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%
%     featx=rEEG(x,Fs,'rEEG_width');
%     
%     
%
% [1] D O’Reilly, MA Navakatikyan, M Filip, D Greene, & LJ Van Marter (2012). Peak-to-peak
% amplitude in neonatal brain monitoring of premature infants. Clinical Neurophysiology,
% 123(11), 2139–53.

% John M. O' Toole, University College Cork
% Started: 19-04-2016
%
% last update: Time-stamp: <2017-03-14 16:21:44 (otoolej)>
%-------------------------------------------------------------------------------
function [featx,reeg_all]=rEEG(x,Fs,feat_name,params_st)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='rEEG_mean'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end
reeg_all=[];


if(isempty(params_st))
    neural_parameters;
    if(strfind(feat_name,'rEEG'))
        params_st=feat_params_st.rEEG;
    else
        params_st=feat_params_st.(char(feat_name));
    end
end

freq_bands=params_st.freq_bands;
N_freq_bands=size(freq_bands,1);
if(isempty(freq_bands))
    N_freq_bands=1;
end



for n=1:N_freq_bands

    % filter (if necessary)
    if(~isempty(freq_bands))
        x_filt=filter_butterworth_withnans(x,Fs,freq_bands(n,2),freq_bands(n,1),5,...
                                           params_st.FILTER_REPLACE_ARTEFACTS);
    else        
        x_filt=x;
    end

    % generate rEEG
    reeg=gen_rEEG(x_filt,Fs,params_st.L_window,params_st.overlap, ...
                  params_st.window_type,params_st.APPLY_LOG_LINEAR_SCALE);
    N=length(reeg);

    DBplot=0;
    if(DBplot)
        figure(1); clf; hold all;
        ttime=0:(length(reeg)-1); ttime=ttime.*params_st.L_window;
        plot_rEEG(ttime,reeg);
        grid on;
    end

    
    switch feat_name
      case 'rEEG_mean'
        %---------------------------------------------------------------------
        % mean rEEG
        %---------------------------------------------------------------------
        featx(n)=nanmean(reeg);
        
      case 'rEEG_median'
        %---------------------------------------------------------------------
        % mean rEEG
        %---------------------------------------------------------------------
        featx(n)=nanmedian(reeg);

      case 'rEEG_lower_margin'
        %---------------------------------------------------------------------
        % 5th prcentile rEEG
        %---------------------------------------------------------------------
        featx(n)=prctile(reeg,5);

      case 'rEEG_upper_margin'
        %---------------------------------------------------------------------
        % 95th prcentile rEEG
        %---------------------------------------------------------------------
        featx(n)=prctile(reeg,95);

      case 'rEEG_width'
        %---------------------------------------------------------------------
        % amplitude bandwidth 
        %---------------------------------------------------------------------
        featx(n)=prctile(reeg,95) - prctile(reeg,5);

      case 'rEEG_SD'
        %---------------------------------------------------------------------
        % standard deviation
        %---------------------------------------------------------------------
        featx(n)=nanstd(reeg);

      case 'rEEG_CV'
        %---------------------------------------------------------------------
        % coefficient of variation 
        %---------------------------------------------------------------------
        featx(n)=nanstd(reeg)/nanmean(reeg);

      case 'rEEG_asymmetry'
        %---------------------------------------------------------------------
        % coefficient of variation 
        %---------------------------------------------------------------------
        A=nanmedian(reeg) - prctile(reeg,5);
        B=prctile(reeg,95) - nanmedian(reeg);
        featx(n)=(B-A)/(A+B);
        
      otherwise
        fprintf('unknown feature ''%s''; check spelling\n',feat_name);
        featx=NaN;
    end

    if(nargout>1)
        reeg_all(n,:)=reeg;
    end
end

    
    
function reeg=gen_rEEG(x,Fs,win_length,win_overlap,win_type,APPLY_LOG_LINEAR_SCALE)
%---------------------------------------------------------------------
% generate the peak-to-peak measure (rEEG)
%---------------------------------------------------------------------
[L_hop,L_epoch,win_epoch]=gen_epoch_window(win_overlap,win_length,win_type,Fs);


N=length(x);
N_epochs=floor( (N-(L_epoch-L_hop))/L_hop );
if(N_epochs<1) N_epochs=1; end
nw=0:L_epoch-1;

%---------------------------------------------------------------------
% generate short-time FT on all data:
%---------------------------------------------------------------------
reeg=NaN(1,N_epochs);
for k=1:N_epochs
    nf=mod(nw+(k-1)*L_hop,N);
    x_epoch=x(nf+1).*win_epoch(:)';

    reeg(k)=max(x_epoch)-min(x_epoch);
end
% no need to resample (as per [1])


% log--linear scale:
if(APPLY_LOG_LINEAR_SCALE)
    ihigh=find(reeg>50);
    if(~isempty(ihigh))
        reeg(ihigh)=50.*log(reeg(ihigh))./log(50);
    end
end


function plot_rEEG(t1,reeg)
%---------------------------------------------------------------------
% plot the linear -log scaled rEEG
%---------------------------------------------------------------------
col_reeg={[0.2539 0.4102 0.8789]};


% log-linear scale:
ihigh=find(reeg>50);
if(~isempty(ihigh))
    reeg(ihigh)=50.*log(reeg(ihigh))./log(50);
end

plot(t1,reeg,'color',col_reeg{1});
xlim([-0.1 t1(end)+0.1]);


ytics=[0 25 50 50*log(250)/log(50) 50*log(1000)/log(50)];
log_ytics_labels={'0','25','50','250','1000'};
set(gca,'ytick',ytics,'yticklabel',log_ytics_labels);
ylim([ytics(1) ytics(end)+5]);


for n=1:length(ytics)
    lp=line(xlim,[ytics(n) ytics(n)],'color',[1 1 1].*0.5);
    uistack(lp,'bottom');
end
