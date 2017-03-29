%-------------------------------------------------------------------------------
% connectivity_features: connectivity features
%
% Syntax: featx=connectivity_features(x,Fs,feat_name,params_st)
%
% Inputs: 
%     x          - epoch of EEG data (size 1 x N)
%     Fs         - sampling frequency (in Hz)
%     feat_name  - feature type, defaults to 'connectivity_BSI';
%                  see full list of 'connectivity_' features in all_features_list.m
%     params_st  - parameters (as structure); 
%                  see neural_parameters.m for examples
%
% Outputs: 
%     featx  - feature at each frequency band 
%
% Example:
%     Fs=64; 
%     data_st=gen_test_EEGdata(32,Fs,1);
%     x=data_st.eeg_data;
%     channel_labels=data_st.ch_labels;
%
%     featx=connectivity_features(x,Fs,'connectivity_corr',[],channel_labels);
%     
%
% [1] van Putten, MJAM (2007). The revised brain symmetry index. Clinical Neurophysiology,
%     118(11), 2362–2367.
% [2] Prichard D, Theiler J (1994). Generating surrogate data for time series with several
%     simultaneously measured variables. Physical Review Letters, 1994;73(7):951–954.
% [3] Faes L, Pinna GD, Porta A, Maestri R, Nollo G (2004). Surrogate data analysis for
%     assessing the significance of the coherence function. IEEE Transactions on
%     Biomedical Engineering, 51(7):1156–1166.
%

% John M. O' Toole, University College Cork
% Started: 13-04-2016
%
% last update: Time-stamp: <2017-03-29 09:26:29 (otoolej)>
%-------------------------------------------------------------------------------
function featx=connectivity_features(x,Fs,feat_name,params_st,ch_labels)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='connectivity_BSI'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end
if(nargin<5 || isempty(ch_labels)), ch_labels=[]; end


DBplot=1;


[N_channels,N]=size(x);
if(N_channels<2), error('require at least 2 channels'); end


% so far no parameters but maybe later:
if(isempty(params_st))
    neural_parameters;
    if(strfind(feat_name,'connectivity'))
        params_st=feat_params_st.connectivity;
    else
        params_st=feat_params_st.(char(feat_name));
    end
end
freq_bands=params_st.freq_bands;

N_freq_bands=size(freq_bands,1);
if(isempty(freq_bands))
    N_freq_bands=1;
end


if(N_channels>2  && ~isempty(ch_labels))
    [ileft,iright]=channel_hemispheres(ch_labels);
    ipairs=channel_hemisphere_pairs(ch_labels);    
elseif(N_channels==2)
    % if no channel labels then guess:
    ileft=1; iright=2;
    ipairs=[1 2]';
end

    
switch feat_name
  case 'connectivity_BSI'
    %---------------------------------------------------------------------
    % brain sysmetry index (revised version, Van Putten, 2007)
    %
    % van Putten, M. J. A. M. (2007). The revised brain symmetry index. Clinical
    % Neurophysiology, 118(11), 2362–2367. http://doi.org/10.1016/j.clinph.2007.07.019
    %---------------------------------------------------------------------

    % a) PSD estimate (Welch's periodogram):
    for k=1:N_channels
        x_epoch=x(k,:);
        x_epoch(isnan(x_epoch))=[];
        
        params_st.method='PSD';
        X(k,:)=gen_spectrum(x_epoch,Fs,params_st,1);
    end

    N=size(X,2); Nfreq=2*(N-1); f_scale=(Nfreq/Fs);

    if(length(ileft)>1)
        X_left=nanmean(X(ipairs(1,:),:),1);
        X_right=nanmean(X(ipairs(2,:),:),1);            
    else
        X_left=X(ileft,:);
        X_right=X(iright,:);            
    end
    
        
    if(DBplot)
        fpp=linspace(0,Fs/2,N);
        fp=fpp;
        set_figure(1); clf; hold all;
        subplot(321); hold all;
        plot(fpp,10*log10(X(ileft,:)'))
        subplot(322); hold all;
        plot(fpp,10*log10(X(iright,:)'))
        subplot(3,2,3:4); hold all;
        plot(fpp,10*log10(X_left'));
        plot(fpp,10*log10(X_right'));
        subplot(3,2,5:6); hold all;
        plot(fpp,(abs( (X_left - X_right)./(X_left + X_right) )));
        disp('--- paused; hit key to continue ---'); pause;
    end

    
    for n=1:N_freq_bands
        ibandpass=ceil(freq_bands(n,1)*f_scale):floor(freq_bands(n,2)*f_scale);        
        ibandpass=ibandpass+1;
        ibandpass(ibandpass<1)=1; ibandpass(ibandpass>N)=N;    
        

        featx(n)=nanmean(abs( (X_left(ibandpass) - X_right(ibandpass)) ./ ...
                              (X_left(ibandpass) + X_right(ibandpass)) ));

        if(DBplot)
            line([1 1].*fp(ibandpass(1)),ylim,'color','k');
            line([1 1].*fp(ibandpass(end)),ylim,'color','k');        
        end
    end

    
  case {'connectivity_coh_mean','connectivity_coh_max','connectivity_coh_freqmax'}
    %---------------------------------------------------------------------
    % coherence (using Welch's PSD)
    %---------------------------------------------------------------------

    % PSD and x-PSD parameters
    N_pairs=size(ipairs,2);

    DBplot=1;    
    if(DBplot), figure(1); clf; hold all;  end
    
    
    featx_pairs=NaN(N_freq_bands,N_pairs);
    for p=1:N_pairs        
        x1=x(ipairs(1,p),:);
        x2=x(ipairs(2,p),:);
        x1(isnan(x1))=[];
        x2(isnan(x2))=[];        
        
        [coh(p,:),pxx,pyy,pxy]=gen_coherence(x1,x2,Fs,params_st);        
        
        set_figure(3); clf; hold all;
        plot(10*log10(pxx)); plot(10*log10(pyy)); plot(10*log10(coh(p,:)));
        plot(plot(log10(abs(pxy).^2)));
        keyboard;
        disp('--- paused; hit key to continue ---'); pause;

        
        % if generating a null-hypothesis distribution from surrogate data:
        if(params_st.coherence_surr_data>0)
            L_surr=params_st.coherence_surr_data;
            coh_surr=zeros(L_surr,size(coh,2));
            
            % generate surrogate signals:
            x1_surr=rand_phase(x1,L_surr);
            x2_surr=rand_phase(x2,L_surr);                

            for m=1:L_surr
                coh_surr(m,:)=gen_coherence(x1_surr(m,:),x2_surr(m,:),Fs,params_st, ...
                                            pxx,pyy);
            end

            % because are not generating pxx and pyy for every iteration (to save
            % compution time), then coherence function may not be bounded to [0,1]:
            coh_surr(coh_surr<0)=0;  coh_surr(coh_surr>1)=1;
            
            % estimating frequency-dependent threshold at the p<α level of significance
            coh_thres=prctile(coh_surr,100*(1-params_st.coherence_surr_alpha));
        end
        dispVars(L_surr,median(coh_thres));

        % if using surrogate data, zero anything below the threhold:
        if(params_st.coherence_surr_data>0)
            coh(p,coh(p,:)<coh_thres)=0;
        end
        
        
        
        if(DBplot)
            set_figure(1);
            subplot(2,2,p); hold all; 
            fp=linspace(0,Fs/2,length(coh(p,:)));
            plot(fp,coh(p,:));
            if(params_st.coherence_surr_data>0)
                plot(fp,coh_thres,'linewidth',2,'color','k');
            end
            xlim([0 fp(end)]);
        end
        
        
        N=size(coh,2); Nfreq=2*(N-1); f_scale=(Nfreq/Fs);
        fp=linspace(0,Fs/2,N);        
        for n=1:N_freq_bands
            ibandpass=ceil(freq_bands(n,1)*f_scale):floor(freq_bands(n,2)*f_scale);        
            ibandpass=ibandpass+1;
            ibandpass(ibandpass<1)=1; ibandpass(ibandpass>N)=N;    
            
            if(strcmp(feat_name,'connectivity_coh_mean'))
                featx_pairs(n,p)=nanmean(coh(p,ibandpass));
            elseif(strcmp(feat_name,'connectivity_coh_max'))
                featx_pairs(n,p)=max(coh(p,ibandpass));
            elseif(strcmp(feat_name,'connectivity_coh_freqmax'))
                [~,imax]=max(coh(p,ibandpass));
                featx_pairs(n,p)=fp(ibandpass(imax));
            end
            
            if(DBplot)
                line([1 1].*fp(ibandpass(1)),ylim,'color','k');
                line([1 1].*fp(ibandpass(end)),ylim,'color','k');        
            end
        end
    end
    
    % median value across channel pairs:
    featx=nanmedian(featx_pairs,2);
    
    
    
  case 'connectivity_corr'
    %---------------------------------------------------------------------
    % cross-correlation (Pearson)
    %---------------------------------------------------------------------
    N_pairs=size(ipairs,2);
    for n=1:N_freq_bands
        for p=1:N_channels
            [x_filt(p,:),inans{p}]=filter_butterworth_withnans(x(p,:),Fs,freq_bands(n,2), ...
                                                    freq_bands(n,1),5, ...
                                                    params_st.FILTER_REPLACE_ARTEFACTS);
        end

        cc_pairs=NaN(1,N_pairs);
        for p=1:N_pairs
            all_inans=unique([inans{ipairs(1,p)} inans{ipairs(2,p)}]);
            
            x1=x_filt(ipairs(1,p),:);
            x2=x_filt(ipairs(2,p),:);
            if(~isempty(all_inans))
        	x1(all_inans)=[];
        	x2(all_inans)=[];                
            end
            
            env1=abs( hilbert(x1) ).^2;
            env2=abs( hilbert(x2) ).^2;                
            
            cc_pairs(p)=corr(env1',env2','type','pearson');
        end
        featx(n)=nanmedian(cc_pairs);
    end

    

  case 'connectivity_lag_corr'
    %---------------------------------------------------------------------
    % lag of cross-correlation 
    % (NOT FINISHED!! needs more testing)
    %---------------------------------------------------------------------
% $$$     N_pairs=size(ipairs,2);    
% $$$     for n=1:N_freq_bands
% $$$         for p=1:N_channels
% $$$             [x_filt(p,:),inans{p}]=filter_butterworth_withnans(x(p,:),Fs,freq_bands(n,2), ...
% $$$                                                               freq_bands(n,1),5, ...
% $$$                                                               params_st.FILTER_REPLACE_ARTEFACTS);
% $$$         end
% $$$ 
% $$$         time_max_lag=NaN(1,N_pairs);
% $$$         for p=1:N_pairs
% $$$             all_inans=unique([inans{ipairs(1,p)} inans{ipairs(2,p)}]);
% $$$             
% $$$             x1=x_filt(ipairs(1,p),:);
% $$$             x2=x_filt(ipairs(2,p),:);
% $$$             if(~isempty(all_inans))
% $$$         	x1(all_inans)=[];
% $$$         	x2(all_inans)=[];                
% $$$             end
% $$$         
% $$$         
% $$$             [cc,lag]=xcorr(x1,x2,'biased');
% $$$ 
% $$$             [~,imax]=max(abs(cc));
% $$$             time_max_lag(p)=lag(imax)/Fs;
% $$$         end
% $$$ 
% $$$         featx(n)=nanmean(time_max_lag);
% $$$ 
% $$$     end
    
    
  case 'connectivity_asynchrony'
    %---------------------------------------------------------------------
    % inter-hemispheric asynchrony (Korotchikova et al. 2011)
    % NOT FINISHED: need to update BSI estimate (see above)
    %---------------------------------------------------------------------
% $$$     for n=1:N_freq_bands
% $$$ 
% $$$         
% $$$         % A) BSI:
% $$$         X_left=abs(fft(x(1,:))).^2; X_right=abs(fft(x(2,:))).^2;
% $$$         X_left=X_left(ibandpass);   X_right=X_right(ibandpass);        
% $$$ 
% $$$         bsi=1-nanmean(abs( (X_left - X_right)./(X_left + X_right) ));
% $$$ 
% $$$         % B) correlation coefficient:
% $$$         x_filt(1,:)=filt_butterworth(x(1,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
% $$$         x_filt(2,:)=filt_butterworth(x(2,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
% $$$         
% $$$         cc=abs(corr(x_filt(1,:)',x_filt(2,:)'));
% $$$ 
% $$$ % $$$         dispVars(cc,bsi);
% $$$         featx(n)=(cc-bsi);
% $$$         
% $$$     end
    
  otherwise
    fprintf('unknown feature ''%s''; check spelling\n',feat_name);
    featx=NaN;
end



function [c,pxx,pyy,pxy]=gen_coherence(x,y,Fs,param_st,pxx,pyy)
%---------------------------------------------------------------------
% generate coherence (magnitude only) between x and y
%---------------------------------------------------------------------
if(nargin<5 || isempty(pxx)), pxx=[]; end

if(isempty(pxx))
    pxx=gen_spectrum(x,Fs,param_st);
    pyy=gen_spectrum(y,Fs,param_st);
end
pxy=gen_cross_spectrum(x,y,Fs,param_st);

c=(abs(pxy).^2)./(pxx.*pyy);


