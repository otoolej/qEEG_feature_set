%-------------------------------------------------------------------------------
% remove_artefacts: simple procedure to remove artefacts
%
% Syntax: data=remove_artefacts(data,ch_labels,Fs,data_ref,ch_refs)
%
% Inputs: 
%     data      - EEG data, in bipolar montage; size: N_channels x N
%     ch_labels - cell of bipolar channel labels, 
%                 e.g. {'C3-O1','C4-O2', 'F3-C3', 'F4-C4'}
%     Fs        - sampling frequency (in Hz)
%     data_ref  - EEG data, in referential  montage; size: (N_channels+1) x N
%     ch_refs   - cell of referential channel labels, 
%                 e.g. {'C3','C4','F3','F4'}
%
% Outputs: 
%     data - EEG data after processing, in bipolar montage, size: N_channels x N
%
% Example:
%     Fs=256;
%     data_st=gen_test_EEGdata(2*60,Fs,1);
%     N=size(data_st.eeg_data_ref,2);
% 
%     % simulate artefacts:
%     % 1. F3 not properly attached:
%     if3=find(strcmp(data_st.ch_labels_ref,'F3'));
%     data_st.eeg_data_ref(if3,:)=randn(1,N).*10;
%
%     % 2. electrode coupling between C4 and Cz
%     ic4=find(strcmp(data_st.ch_labels_ref,'C4'));
%     icz=find(strcmp(data_st.ch_labels_ref,'Cz'));
%     data_st.eeg_data_ref(icz,:)=data_st.eeg_data_ref(ic4,:)+randn(1,N).*5;
%
%     % re-generate bipolar montage:
%     [data_st.eeg_data,data_st.ch_labels] = ...
%        set_bi_montage(data_st.eeg_data_ref,data_st.ch_labels_ref, ...
%                                 data_st.ch_labels_bi);
%
%     % remove channels:
%     eeg_art=remove_artefacts(data_st.eeg_data,data_st.ch_labels,data_st.Fs, ...
%                         data_st.eeg_data_ref,data_st.ch_labels_ref);


% John M. O' Toole, University College Cork
% Started: 05-04-2016
%
% last update: Time-stamp: <2019-09-05 13:59:26 (otoolej)>
%-------------------------------------------------------------------------------
function data=remove_artefacts(data,ch_labels,Fs,data_ref,ch_refs)
if(nargin<3), error('requires 3 input arguments.'); end
if(nargin<4 || isempty(data_ref)), data_ref=[]; end
if(nargin<5 || isempty(ch_refs)), ch_refs=[]; end


neural_parameters;

DBverbose=0;

[N_channels,N]=size(data);

%---------------------------------------------------------------------
% 0. check in referential mode first; is there problem with one 
%    channel (e.g. Cz)
%---------------------------------------------------------------------
irem_channel=[];
if(~isempty(data_ref))
    x_filt=zeros(size(data_ref));
    for n=1:size(data_ref,1)
        x_filt(n,:)=filter_butterworth_withnans(data_ref(n,:),Fs,20,0.5,5);
    end

    
    r=corrcoef(x_filt');
    r(1:size(r,1)+1:end)=NaN;
    r_channel=nanmean(r);
    clear x_filt;

    ilow=find(abs(r_channel)<ART_REF_LOW_CORR);
    if(~isempty(ilow))
        nn=1; irem_channel=[];
        for n=ilow
            ch_find=upper(ch_refs{n});
            itmp=find(cellfun(@(x) ~isempty(strfind(x,ch_find)), upper(ch_labels)));
            
            irem_channel=[irem_channel itmp];
            nn=nn+1;
        end
        
        fprintf(':: remove channel (low ref. correlation): %s\n',ch_labels{irem_channel});
    end
    if(DBverbose)
        print_table(r_channel',{'corr'},ch_refs)    
    end
    data(irem_channel,:)=NaN;
end

ichannels=1:N_channels;
ichannels(irem_channel)=[];
N_channels=length(ichannels);


%---------------------------------------------------------------------
% 1. look for electrode coupling:
%---------------------------------------------------------------------
A=[];
if(N_channels>4)
    [ileft,iright]=channel_hemispheres(ch_labels(ichannels));

    if(length(ileft)>1 && length(iright)>1)
        x_means=zeros(1,N_channels);

        x_filt=zeros(N_channels,N);
        for n=1:N_channels
            x_filt(n,:)=filter_butterworth_withnans(data(ichannels(n),:),Fs,20,0.5,5);
        end


        for n=1:N_channels
            x_means(n)=nanmean( abs( x_filt(n,:) ).^2 );
            
            A{n,1}=x_means(n);
            A{n,2}=ch_labels{ichannels(n)};
        end
        clear x_filt;

        % 1/4 of the median of channel energy:
        cut_off_left=median(x_means(ileft))/4;
        cut_off_right=median(x_means(iright))/4;


        ishort_left=find(x_means(ileft)<cut_off_left);
        ishort_right=find(x_means(iright)<cut_off_right);

        ishort=[ileft(ishort_left) iright(ishort_right)];
        ishort=ichannels(ishort);
        if(~isempty(ishort))
            fprintf(':: remove channel (electrode coupling): %s\n',ch_labels{ishort});
            data(ishort,:)=NaN;
            irem_channel=[irem_channel ishort];
        end
    end
end

if(DBverbose && ~isempty(A)),  print_table(A); end



ichannels=1:size(data,1);
ichannels(irem_channel)=[];
N_channels=length(ichannels);



% all other artefacts are on a channel-by-channel basis:
irem=[];
for n=ichannels
    data(n,:)=art_per_channel(data(n,:),Fs);
    
    irem=[irem find(isnan(data(n,:)))];
end

% remove artefacts across all channels
data(ichannels,unique(irem))=NaN;



function x=art_per_channel(x,Fs)
%---------------------------------------------------------------------
% remove artefacts on a per-channel basis
%---------------------------------------------------------------------
neural_parameters;

DBverbose=1;

N=length(x);

%---------------------------------------------------------------------
% 1. electrode-checks (continuous row of zeros)
%---------------------------------------------------------------------
x_channel=x;
x_channel(x_channel~=0)=1;
irem=zeros(1,N);    
[lens,istart,iend]=len_cont_zeros(x_channel,0);
ielec=find(lens>=(ART_ELEC_CHECK*Fs));
if(~isempty(ielec) && DBverbose)
    fprintf(['electrode check at time(s): ' repmat('%g ',1,length(ielec))  ...
             ' \n'],istart(ielec)./Fs);
end
for m=ielec
    irun=[istart(m)-1:iend(m)+1];
    irun(irun<1)=1; irun(irun>=N)=N;
    irem(irun)=1;
    x(irun)=NaN;
end
if(any(irem==1) && DBverbose)
    fprintf('continuous row of zeros: %.2f%%\n', ...
            100*length(find(irem==1))/length(x));
end

x_nofilt=x;
[x_filt,inans]=filter_butterworth_withnans(x,Fs,40,0.1,[5 2]);


%---------------------------------------------------------------------
% 2. high-amplitude artefacts
%---------------------------------------------------------------------
art_coll=ART_TIME_COLLAR*Fs;
irem=zeros(1,N);    

x_hilbert=abs( hilbert(x_filt) );    

thres_upper=ART_HIGH_VOLT;
ihigh=find(x_hilbert>thres_upper);
if(~isempty(ihigh))
    for p=1:length(ihigh)
        irun=(ihigh(p)-art_coll):(ihigh(p)+art_coll);
        irun(irun<1)=1;  irun(irun>N)=N;               
        irem(irun)=1;
    end
end
x(irem==1)=NaN;
if(any(irem==1) && DBverbose)
    fprintf('length of high-amplitude artefacts: %.2f%%\n', ...
            100*length(find(irem==1))/length(x));
end




%---------------------------------------------------------------------
% 3. continuous constant values (i.e. artefacts)
%---------------------------------------------------------------------
art_coll=ART_DIFF_TIME_COLLAR*Fs;
x_diff_all=zeros(1,N);
irem=zeros(1,N);    

x_diff_all=[diff(x) 0];        
x_diff=x_diff_all;
x_diff(x_diff~=0)=1;
[lens,istart,iend]=len_cont_zeros(x_diff,0);

% if exactly constant for longer than . then remove:
ielec=find(lens>=(ART_DIFF_MIN_TIME*Fs));
% $$$ if(~isempty(ielec))
% $$$     fprintf(['channel: %s; flat voltage at time(s): ' repmat('%g ',1,length(ielec))  ...
% $$$              ' \n'],bi_mont{n},istart(ielec)./Fs);
% $$$ end
for m=ielec
    irun=[(istart(m)-art_coll):(iend(m)+art_coll)];
    irun(irun<1)=1; irun(irun>=N)=N;
    irem(irun)=1;
    x(irun)=NaN;
end
if(any(irem==1) && DBverbose)
    fprintf('continuous row of constant values: %.2f%%\n', ...
            100*length(find(irem==1))/length(x));
end


%---------------------------------------------------------------------
% 4. sudden jumps in amplitudes or constant values (i.e. artefacts)
%---------------------------------------------------------------------
art_coll=ART_DIFF_TIME_COLLAR*Fs;
irem=zeros(1,N);    
x_diff=x_diff_all;    
    
ihigh=find(abs(x_diff)>ART_DIFF_VOLT);
if(~isempty(ihigh))
    for p=1:length(ihigh)
        irun=(ihigh(p)-art_coll):(ihigh(p)+art_coll);
        irun(irun<1)=1;  irun(irun>N)=N;               
        irem(irun)=1;
    end
end
xb=x;
x(irem==1)=NaN;


% before filtering, but should be eliminated anyway
x(inans)=NaN;
inans=find(isnan(x));
x_nofilt(inans)=NaN;
x=x_nofilt;


if(any(irem==1) && DBverbose)
    fprintf('length of sudden-jump artefacts: %.2f%%\n', ...
            100*length(find(irem==1))/length(x));
    
% $$$     figure(22); clf; hold all;
% $$$     plot(x); plot(xb);
% $$$     disp('--- paused; hit key to continue ---'); pause;

end

