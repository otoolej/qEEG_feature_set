%---------------------------------------------------------------------
% generate EEG-like data (coloured Gaussian noise)
%---------------------------------------------------------------------
Fs=256;
data_st=gen_test_EEGdata(2*60,Fs,1);


%---------------------------------------------------------------------
% generate artifical artefacts:
%---------------------------------------------------------------------
N=size(data_st.eeg_data_ref,2);

% 1. F4 not properly attached:
if3=find(strcmp(data_st.ch_labels_ref,'F3'));
data_st.eeg_data_ref(if3,:)=randn(1,N).*10;

% 2. electrode coupling between C4 and Cz
ic4=find(strcmp(data_st.ch_labels_ref,'C4'));
icz=find(strcmp(data_st.ch_labels_ref,'Cz'));
data_st.eeg_data_ref(icz,:)=data_st.eeg_data_ref(ic4,:)+randn(1,N).*5;


%---------------------------------------------------------------------
% DO PLOTS:
% only if EEG_viewer installed (https://github.com/otoolej/eeg_viewer.git)
%---------------------------------------------------------------------
do_plot=0;
if(do_plot)
    eeg_plotgui_withannos('signals',data_st.eeg_data_ref, ...
                          'Fs',Fs, ...
                          'channel_labels',data_st.ch_labels_ref, ...
                          'epoch_length',20, ...
                          'insert_ta_scale',1,'ta_xlength',2);
end
%---------------------------------------------------------------------
% .... END PLOTS
%---------------------------------------------------------------------


% re-generate bipolar montage:
[data_st.eeg_data,data_st.ch_labels] = ...
    set_bi_montage(data_st.eeg_data_ref,data_st.ch_labels_ref, ...
                                 data_st.ch_labels_bi);


% remove channels:
eeg_art=remove_artefacts(data_st.eeg_data,data_st.ch_labels,data_st.Fs, ...
                         data_st.eeg_data_ref,data_st.ch_labels_ref);



if(do_plot)
    eeg_plotgui_withannos('signals',eeg_art, ...
                          'Fs',Fs, ...
                          'channel_labels',data_st.ch_labels, ...
                          'bipolar_montage',-1, ...
                          'epoch_length',30);
end





