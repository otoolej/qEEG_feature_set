%-------------------------------------------------------------------------------
% gen_test_EEGdata: generate EEG-like data as coloured Gaussian noise 
%
% Syntax: data_st=gen_test_EEGdata(dur,Fs,include_bipolar,discont_activity)
%
% Inputs: 
%     dur:              duration of EEG-like data in seconds (default 300 seconds)
%     Fs:               sampling frequency
%     include_bipolar:  include the bipolar montage aswell as the referential  
%     discont_activity: discontinuous-like activity of preterm EEG (bursts and inter-bursts)
%
% Outputs: 
%     data_st:          structure including EEG data (referential montage), sampling
%                       frequency, and channel labels
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data(1,:);
%     
%     figure(1); clf; hold all;
%     ttime=(0:(length(x)-1))./Fs;
%     plot(ttime,x);
%     xlabel('time (seconds)'); ylabel('\muV');
%

% John M. O' Toole, University College Cork
% Started: 01-09-2016
%
% last update: Time-stamp: <2017-03-14 17:35:14 (otoolej)>
%-------------------------------------------------------------------------------
function data_st=gen_test_EEGdata(dur,Fs,include_bipolar,discont_activity)
if(nargin<1 || isempty(dur)), dur=60*2; end
if(nargin<2 || isempty(Fs)), Fs=256; end
if(nargin<3 || isempty(include_bipolar)), include_bipolar=0; end
if(nargin<4 || isempty(discont_activity)), discont_activity=0; end


N=floor(dur*Fs);


%---------------------------------------------------------------------
% generate white Gaussian noise and filter with moving average filter
% (referential montage has common signal)
%---------------------------------------------------------------------
L_ma=Fs/4;
L_ma_noise=ceil(Fs/16);
N_channels=9;
eeg_common=5.*filter(ones(1,L_ma)./L_ma,1,randn(1,N)); 

eeg_data_ref=filter(ones(1,L_ma)./L_ma,1,randn(N_channels,N)')' + ...
    (filter(ones(1,L_ma_noise)./L_ma_noise,1,randn(N_channels,N)')')./100 + ...
    repmat(eeg_common,[N_channels 1]);

if(discont_activity)
    % if want discontinuous-like activity (burst and inter-bursts):
    ibursts=randi([1 N],1,fix(dur/100)*Fs);
    amps=abs(1+rand(1,fix(dur/100)*Fs));
    
    bursts=zeros(1,N);
    bursts(ibursts)=1.*amps;
    
    L_win=Fs*2;
    bursts_smooth=100.*filter(hamming(L_win)./L_win,1,bursts); 
    
    eeg_data_ref=10.*bsxfun(@plus,eeg_data_ref,bsxfun(@times,eeg_data_ref, 20.* ...
                                               bursts_smooth));
else
    % slow-duration amplitude modulation:
    % not included as needs more work:
% $$$     long_env=randn(1,N);
% $$$     long_env=long_env-mean(long_env);
% $$$     long_env=filter(ones(1,100*Fs)./(100*Fs),1,long_env);
% $$$     long_env=abs(long_env);
% $$$     long_env=sqrt(.5).*long_env./std(long_env);
% $$$     long_env=(long_env)+0.2;
% $$$     figure(2); clf; hold all;
% $$$     plot(long_env);
% $$$     eeg_data_ref=bsxfun(@times,eeg_data_ref,long_env);

    
    eeg_data_ref=eeg_data_ref.*50;
end


data_st.eeg_data_ref=eeg_data_ref;
data_st.Fs=Fs;


%---------------------------------------------------------------------
% generate bipolar and monopolar labels:
%---------------------------------------------------------------------
data_st.ch_labels_bi={{'F4','C4'},{'F3','C3'},{'C4','T4'},{'C3','T3'}, ...
                   {'C4','Cz'},{'Cz','C3'}, ...
                   {'C4','O2'},{'C3','O1'}};
data_st.ch_labels_ref=unique([data_st.ch_labels_bi{:}]);



%---------------------------------------------------------------------
% include bipolar montage as well?
%---------------------------------------------------------------------
if(include_bipolar)
    % generate bi-polar montage:
    [data_st.eeg_data,data_st.ch_labels]=set_bi_montage(data_st.eeg_data_ref, ...
                                                      data_st.ch_labels_ref, ...
                                                      data_st.ch_labels_bi);
end


% plot:
DBplot=0;
if(DBplot)
    
    %---------------------------------------------------------------------
    % EEG viewer installed?
    % (https://github.com/otoolej/eeg_viewer)
    %---------------------------------------------------------------------
    eeg_plotgui_withannos('signals',data_st.eeg_data_ref, ...
                          'Fs',data_st.Fs, ...
                          'channel_labels',data_st.ch_labels_ref, ...
                          'epoch_length',30);
end



