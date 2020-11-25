%-------------------------------------------------------------------------------
% connectivity_features: connectivity features
%
% Syntax: featx = connectivity_features(x, Fs, feat_name, params_st)
%
% Inputs: 
%     x          - epoch of EEG data (size N_channels x N)
%     Fs         - sampling frequency (in Hz)
%     feat_name  - feature type, defaults to 'connectivity_BSI';
%                  see full list of 'connectivity_' features in all_features_list.m
%     params_st  - parameters (as structure); 
%                  see neural_parameters.m for examples
%     ch_labels  - cell function with channel names (1 x N_channels)
%
% Outputs: 
%     featx  - feature at each frequency band 
%
% Example:
%     Fs = 64; 
%     data_st = gen_test_EEGdata(32, Fs, 1);
%     x = data_st.eeg_data;
%     channel_labels = data_st.ch_labels;
%
%     featx = connectivity_features(x, Fs, 'connectivity_corr', [], channel_labels);
%     
%
% [1] van Putten, MJAM (2007). The revised brain symmetry index. Clinical Neurophysiology,
%     118(11), 2362–2367.
% [2] Prichard D, Theiler J (1994). Generating surrogate data for time series with several
%     simultaneously measured variables. Physical Review Letters, 1994;73(7):951–954.
% [3] Faes L, Pinna GD, Porta A, Maestri R, Nollo G (2004). Surrogate data analysis for
%     assessing the significance of the coherence function. IEEE Transactions on
%     Biomedical Engineering, 51(7):1156–1166.
% [4] Halliday, DM, Rosenberg, JR, Amjad, AM, Breeze, P, Conway, BA, &
%     Farmer, SF. (1995). A framework for the analysis of mixed time series/point
%     process data--theory and application to the study of physiological tremor, single
%     motor unit discharges and electromyograms. Progress in Biophysics and Molecular
%     Biology, 64(2–3), 237–278.
%

% John M. O' Toole, University College Cork
% Started: 13-04-2016
%
% last update: Time-stamp: <2020-11-25 19:14:45 (otoolej)>
%-------------------------------------------------------------------------------
function featx = connectivity_features(x, Fs, feat_name, params_st, ch_labels)
if(nargin<2), error('need 2 input arguments'); end
if(nargin<3 || isempty(feat_name)), feat_name = 'connectivity_BSI'; end
if(nargin<4 || isempty(params_st)), params_st = []; end
if(nargin<5 || isempty(ch_labels)), ch_labels = []; end


DBplot = 0;
featx = NaN;


[N_channels, N] = size(x);
if(N_channels<2)
    warning('REQUIRED: at least 2 channels for connectivity functions'); 
    return;
end


% so far no parameters but maybe later:
if(isempty(params_st))
    neural_parameters;
    if(strfind(feat_name, 'connectivity'))
        params_st = feat_params_st.connectivity;
    else
        params_st = feat_params_st.(char(feat_name));
    end
end
freq_bands = params_st.freq_bands;

N_freq_bands = size(freq_bands, 1);
if(isempty(freq_bands))
    N_freq_bands = 1;
end

if(~isempty(ch_labels))
    % find left and right channels:
    [ileft, iright] = channel_hemispheres(ch_labels);
    ipairs = channel_hemisphere_pairs(ch_labels);    
else
    warning('REQUIRED: channel names for connectivity functions.');
    return;
end
if(isempty(ipairs))
    warning('no channel pairs (left / right) for connectivity functions.');    
    return;
end




switch feat_name
  case 'connectivity_BSI'
    %---------------------------------------------------------------------
    % brain sysmetry index (revised version, Van Putten, 2007)
    %
    % van Putten, MJAM (2007). The revised brain symmetry index. Clinical
    % Neurophysiology, 118(11), 2362–2367. http://doi.org/10.1016/j.clinph.2007.07.019
    %---------------------------------------------------------------------

    % a) PSD estimate (Welch's periodogram):
    for k = 1:N_channels
        x_epoch = x(k, :);
        x_epoch(isnan(x_epoch)) = [];
        
        [X(k, :), ~, f_scale, ~, fp] = gen_spectrum(x_epoch, Fs, params_st, 1);
    end

    N = size(X, 2); %Nfreq = 2*(N-1); f_scale = (Nfreq/Fs);

    if(length(ileft)>1)
        X_left = nanmean(X(ipairs(1, :), :), 1);
        X_right = nanmean(X(ipairs(2, :), :), 1);            
    else
        X_left = X(ileft, :);
        X_right = X(iright, :);            
    end
    
    
    if(DBplot)
        set_figure(1); clf; hold all;
        subplot(321); hold all;
        plot(fp, 10*log10(X(ileft, :)'))
        subplot(322); hold all;
        plot(fp, 10*log10(X(iright, :)'))
        subplot(3, 2, 3:4); hold all;
        plot(fp, 10*log10(X_left'));
        plot(fp, 10*log10(X_right'));
        subplot(3, 2, 5:6); hold all;
        plot(fp, (abs( (X_left - X_right)./(X_left + X_right) )));
    end

    
    for n = 1:N_freq_bands
        if(n == 1)
            istart = ceil(freq_bands(n, 1)*f_scale);
        else
            istart = ibandpass(end)-1;
        end
        
        ibandpass = istart:floor(freq_bands(n, 2)*f_scale);        
        ibandpass = ibandpass+1;
        ibandpass(ibandpass<1) = 1; ibandpass(ibandpass>N) = N;    
        

        featx(n) = nanmean(abs( (X_left(ibandpass) - X_right(ibandpass)) ./ ...
                                (X_left(ibandpass) + X_right(ibandpass)) ));

        if(DBplot)
            line([1 1].*fp(ibandpass(1)), ylim, 'color', 'k');
            line([1 1].*fp(ibandpass(end)), ylim, 'color', 'k');        
        end
    end

    
  case {'connectivity_coh_mean', 'connectivity_coh_max', 'connectivity_coh_freqmax'}
    %---------------------------------------------------------------------
    % coherence (using Welch's PSD)
    %---------------------------------------------------------------------
    
    % check to see if using the right PSD estimate
    % (only 'PSD' or 'bartlett-PSD' allowed):
    if(strcmp(params_st.method, 'periodogram') || strcmp(params_st.method, 'robust-PSD'))
        fprintf('----------- -WARNING- ------------\n' );
        fprintf('Must use averaging PSD estimate (e.g. Welch PSD) for coherence.\n');
        fprintf('To do so, set: feat_params_st.connectivity.method=''PSD''\n');
        fprintf('in ''neural_parameters.m'' file.\n');
        
        if(strcmp(lower(params_st.coherence_zero_level), 'analytic'))
            params_st.method = 'bartlett-PSD';
            
            warning('Forcing PSD method for connectivity analysis to Bartlett PSD');
            
        else
            warning(sprintf('%s\n%s', ...
                            'Forcing PSD method for connectivity analysis to Welch PSD ', ...
                            '(may need to adjust window type, window size, and overlap.)'));
                   
            params_st.method = 'PSD';
        end
        fprintf('----------------------------------\n' );
    end
    if(strcmp(lower(params_st.coherence_zero_level), 'analytic') && ...
       ~strcmp(params_st.method, 'bartlett-PSD'))
        
        % check if PSD method is Bartlett or not:
        fprintf('----------- -WARNING- ------------\n' );
        fprintf('If want to use the analytic zero-level threshold for\n');
        fprintf('coherence (Halliday et al. 1995) then need to use Barlett PSD.\n');
        fprintf('To do so, set: feat_params_st.connectivity.method=''bartlett-PSD''\n');
        fprintf('in ''neural_parameters.m'' file.\n');
        warning(['Forcing PSD method for connectivity analysis to Bartlett ' ...
                 'PSD']);
        fprintf('----------------------------------\n' );
        
        params_st.method = 'bartlett-PSD';
    end

    
    
    % PSD and x-PSD parameters
    N_pairs = size(ipairs, 2);

    DBplot = 0;    
    if(DBplot), set_figure(1); end
    
    
    featx_pairs = NaN(N_freq_bands, N_pairs);
    for p = 1:N_pairs        
        x1 = x(ipairs(1, p), :);
        x2 = x(ipairs(2, p), :);
        x1(isnan(x1)) = [];
        x2(isnan(x2)) = [];        

        
        % 1) estimate the coherence function:
        [coh(p, :), pxx, pyy, f_scale, fp] = gen_coherence(x1, x2, Fs, params_st); 
        
        
        % 2) if estimating a zero-level threshold for the coherence:
        switch lower(params_st.coherence_zero_level)
          case 'surr'
            %---------------------------------------------------------------------
            % if generating a null-hypothesis distribution from surrogate data:
            %---------------------------------------------------------------------

            L_surr = params_st.coherence_surr_iter;
            coh_surr = zeros(L_surr, size(coh, 2));
            
            % generate surrogate signals:
            x1_surr = rand_phase(x1, L_surr);
            x2_surr = rand_phase(x2, L_surr);                

            for m = 1:L_surr
                coh_surr(m,:) = gen_coherence(x1_surr(m,:), x2_surr(m,:), Fs, params_st); 
            end

            % estimating frequency-dependent threshold at the p<α level of significance
            coh_thres = prctile(coh_surr, 100*(1-params_st.coherence_zero_alpha));
            
          case  'analytic'
            %---------------------------------------------------------------------
            % or if using an analytic method
            %---------------------------------------------------------------------
            % number of segments (no overlap with Bartlett PSD)
            L = floor(length(x1) / (Fs * params_st.L_window));
            
            coh_thres = 1 - (params_st.coherence_zero_alpha) ^ (1 / (L-1) );
            
          case ''
            coh_thres = [];
            
          otherwise
            error(['unknown option for ''coherence_zero_level''; see neural_parameters.m' ...
                   'for details']);
        end
        
        
        % 3) threhold:
        if(~isempty(coh_thres))
            coh(p, coh(p, :)<coh_thres) = 0;
        end
        
        
        if(DBplot)
            % subplot(2, 2, p); hold all; 
            plot(fp, coh(p, :));
            if(~isempty(coh_thres))
                if(length(coh_thres) == 1)
                    line(xlim,  coh_thres .* ones(1, 2), 'linewidth', 2, 'color', 'k');
                else
                    plot(fp, coh_thres, 'linewidth', 2, 'color', 'k');
                end
            end
            xlim([0 fp(end)]);
        end
        
        % compute the coherence (either mean, max, or max. frequency) for each frequency band:
        N = size(coh, 2);
        for n = 1:N_freq_bands
            if(n == 1)
                istart = ceil(freq_bands(n, 1)*f_scale);
            else
                istart = ibandpass(end)-1;
            end
            ibandpass = istart:floor(freq_bands(n, 2)*f_scale);        
            ibandpass = ibandpass+1;
            ibandpass(ibandpass<1) = 1; ibandpass(ibandpass>N) = N;    
            
            if(strcmp(feat_name, 'connectivity_coh_mean'))
                featx_pairs(n, p) = nanmean(coh(p, ibandpass));
            elseif(strcmp(feat_name, 'connectivity_coh_max'))
                featx_pairs(n, p) = max(coh(p, ibandpass));
            elseif(strcmp(feat_name, 'connectivity_coh_freqmax'))
                [~, imax] = max(coh(p, ibandpass));
                featx_pairs(n, p) = fp(ibandpass(imax));
            end
            
            if(DBplot)
                line([1 1].*fp(ibandpass(1)), ylim, 'color', 'k');
                line([1 1].*fp(ibandpass(end)), ylim, 'color', 'k');        
            end
        end
    end
    
    % median value across channel pairs:
    featx = nanmedian(featx_pairs, 2);
    
    
    
  case 'connectivity_corr'
    %---------------------------------------------------------------------
    % cross-correlation (Pearson)
    %---------------------------------------------------------------------
    N_pairs = size(ipairs, 2);
    for n = 1:N_freq_bands
        for p = 1:N_channels
            [x_filt(p, :), inans{p}] = filter_butterworth_withnans(x(p, :), Fs, freq_bands(n, 2),  ...
                                                              freq_bands(n, 1), 5,  ...
                                                              params_st.FILTER_REPLACE_ARTEFACTS);
        end

        cc_pairs = NaN(1, N_pairs);
        for p = 1:N_pairs
            all_inans = unique([inans{ipairs(1, p)} inans{ipairs(2, p)}]);
            
            x1 = x_filt(ipairs(1, p), :);
            x2 = x_filt(ipairs(2, p), :);
            if(~isempty(all_inans))
        	x1(all_inans) = [];
        	x2(all_inans) = [];                
            end
            
            env1 = abs( hilbert(x1) ).^2;
            env2 = abs( hilbert(x2) ).^2;                
            
            cc_pairs(p) = corr(env1', env2', 'type', 'pearson');
        end
        featx(n) = nanmedian(cc_pairs);
    end

    
    
  otherwise
    fprintf('unknown feature ''%s''; check spelling\n', feat_name);
    featx = NaN;
end



function [c, pxx, pyy, f_scale, fp] = gen_coherence(x, y, Fs, param_st)
%---------------------------------------------------------------------
% generate coherence (magnitude only) between x and y
%---------------------------------------------------------------------
pxx = gen_spectrum(x, Fs, param_st);
pyy = gen_spectrum(y, Fs, param_st);
[pxy, ~, f_scale, fp] = gen_cross_spectrum(x, y, Fs, param_st);

c = (abs(pxy).^2)./(pxx.*pyy);


