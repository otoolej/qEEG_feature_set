%-------------------------------------------------------------------------------
% plot_rEEG: plot the linear-log scaled rEEG
%
% Syntax: []=plot_rEEG(t1,reeg)
%
% Inputs: 
%     t1   -  time vector (size 1 x N)
%     reeg -  rEEG (size 1 x N)
%
%
% Example:
%     % first, generate some 2-hours of EEG-like signals:
%     Fs=64; 
%     data_st=gen_test_EEGdata(2*60*60,Fs,1,1);
%
%     % set parameters (or use the defaults in 'neural_parameters.m'):
%     params_st.L_window=2; % in seconds
%     params_st.window_type='rect'; % type of window
%     params_st.overlap=0; % overlap in percentage
%     params_st.APPLY_LOG_LINEAR_SCALE=0; % use this scale (either 0 or 1)
%     params_st.freq_bands=[1 20];
%     params_st.FILTER_REPLACE_ARTEFACTS='nans';
%
%
%     % then, generate the rEEG:
%     [~,reeg]=rEEG(data_st.eeg_data(1,:),Fs,'rEEG_mean',params_st);
%
%
%     % then plot:
%     figure(1); clf; hold all;
%     ttime=0:(length(reeg)-1); ttime=(ttime.*params_st.L_window)/(60*60);
%     plot_rEEG(ttime,reeg);
%     xlabel('time (hours)'); ylabel('voltage');
%

% John M. O' Toole, University College Cork
% Started: 16-06-2017
%
% last update: Time-stamp: <2017-06-16 18:39:47 (otoolej)>
%-------------------------------------------------------------------------------
function plot_rEEG(t1,reeg)
if(nargin<1 || isempty(t1)), t1=[]; end


col_reeg={[0.2539 0.4102 0.8789]};


% log-linear scale:
ihigh=find(reeg>50);
if(~isempty(ihigh))
    reeg(ihigh)=50.*log(reeg(ihigh))./log(50);
end

if(isempty(t1))
    plot(reeg,'color',col_reeg{1});
else
    plot(t1,reeg,'color',col_reeg{1});
    xlim([-0.1 t1(end)+0.1]);
end



ytics=[0 25 50 50*log(250)/log(50) 50*log(1000)/log(50)];
log_ytics_labels={'0','25','50','250','1000'};
set(gca,'ytick',ytics,'yticklabel',log_ytics_labels);
ylim([ytics(1) ytics(end)+5]);


for n=1:length(ytics)
    lp=line(xlim,[ytics(n) ytics(n)],'color',[1 1 1].*0.5);
    uistack(lp,'bottom');
end
