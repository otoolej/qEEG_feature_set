%-------------------------------------------------------------------------------
% connectivity_features: connectivity features
%
% Syntax: featx=connectivity_features(x,Fs,feat_name,params_st)
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

% John M. O' Toole, University College Cork
% Started: 13-04-2016
%
% last update: Time-stamp: <2016-04-26 16:22:47 (otoolej)>
%-------------------------------------------------------------------------------
function featx=connectivity_features(x,Fs,feat_name,params_st,ch_labels)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name='envelope'; end
if(nargin<4 || isempty(params_st)), params_st=[]; end
if(nargin<5 || isempty(ch_labels)), ch_labels=[]; end


DBplot=0;


[N_channels,N]=size(x);
if(N_channels<2), error('require at least 2 channels'); end


% so far no parameters but maybe later:
if(isempty(params_st))
    quant_feats_parameters;
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

x_orig=x;


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
    %---------------------------------------------------------------------

    % a) PSD estimate (Welch's periodogram):
    win_length=make_odd(params_st.PSD_window*Fs);
    overlap=ceil(win_length*(1-params_st.PSD_overlap/100));

    % assuming just two pairs corresponding to left/right:
    X=[]; X_left=[]; X_right=[];
    for k=1:size(x,1)
        [X(k,:),fp]=pwelch(x(k,:),win_length,overlap,[],Fs);
    end

    N=size(X,2); Nfreq=2*(N-1); f_scale=(Nfreq/Fs);

    if(length(ileft)>1)
        X_left=sum(X(ileft,:),1);
        X_right=sum(X(iright,:),1);            
    else
        X_left=X(ileft,:);
        X_right=X(iright,:);            
    end
    
        
    DBplot=0;
    if(DBplot)
        fpp=fp;
        figure(1); clf; hold all;
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
    
    
  case 'connectivity_corr'
    %---------------------------------------------------------------------
    % cross-correlation (Spearmans)
    %---------------------------------------------------------------------
    N_pairs=size(ipairs,2);
    for n=1:N_freq_bands
        for p=1:N_channels
            x_filt(p,:)=filt_butterworth(x(p,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        end

        cc_pairs=NaN(1,N_pairs);
        for p=1:N_pairs
            env1=abs( hilbert(x_filt(ipairs(1,p),:)) ).^2;
            env2=abs( hilbert(x_filt(ipairs(2,p),:)) ).^2;                
            
            cc_pairs(p)=corr(env1',env2','type','pearson');
        end
% $$$         dispVars(cc_pairs);
% $$$         featx(n)=nanmedian(cc_pairs);
    end
    
  case {'connectivity_coh_mean','connectivity_coh_max','connectivity_coh_freqmax'}
    %---------------------------------------------------------------------
    % coherence (using Welch's PSD)
    %---------------------------------------------------------------------

    % PSD and x-PSD parameters
    win_length=make_odd(params_st.PSD_window*Fs);
    overlap=ceil(win_length*(1-params_st.PSD_overlap/100));
    
    N_pairs=size(ipairs,2);

    DBplot=0;    
    if(DBplot)
        figure(1); clf; hold all;    
    end
    
    featx_pairs=NaN(N_freq_bands,N_pairs);
    for p=1:N_pairs        
        pxx=pwelch(x(ipairs(1,p),:),win_length,overlap,[],Fs);
        pyy=pwelch(x(ipairs(2,p),:),win_length,overlap,[],Fs);        
        [pxy,fp]=cpsd(x(ipairs(1,p),:),x(ipairs(2,p),:),win_length,overlap,[],Fs);        
        
        coh(p,:)=(abs(pxy).^2)./(pxx.*pyy);

        if(DBplot)
            subplot(2,2,p); hold all; 
            plot(fp,coh(p,:));
        end
        
        N=size(coh,2); Nfreq=2*(N-1); f_scale=(Nfreq/Fs);
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
    
    

  case 'connectivity_lag_corr'
    %---------------------------------------------------------------------
    % lag of cross-correlation 
    % (not sure about if this one make sense with >2  channels?? )
    %---------------------------------------------------------------------
    for n=1:N_freq_bands
        
        x_filt(1,:)=filt_butterworth(x(1,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        x_filt(2,:)=filt_butterworth(x(2,:),Fs,freq_bands(n,2),freq_bands(n,1),5);
        
        [cc,lag]=xcorr(x_filt(1,:),x_filt(2,:),'biased');

        [~,imax]=max(abs(cc));
        time_max_lag=lag(imax)/Fs;

% $$$         figure(1); clf; hold all;
% $$$         plot(lag,cc);
        
        featx(n)=time_max_lag;

    end
    
    
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
end

