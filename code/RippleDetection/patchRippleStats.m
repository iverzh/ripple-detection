function rippleStats = patchRippleStats(rippleStats, LFPdat, LFPdat_other, sfreq, do_plots, checkSpike)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Applies further rejection criteria to ripples saved in rippleStats
% Notes: 
%   Back up rippleStats before running
%   Run chanSelect before running
%
% INPUTS (all should be post-chanSelect()
%   -LFPdat: broadband cortical LFP, ch x samples
%   -LFPdat_other: broadband LFP data from other sites to apply
%   cross-channel spike rejection, i.e. if LFPdat is NC data then 
%   LFPdat_other might be HC data
%   -rippleStats
%   -do_plots: plot individual examples with lines marking zero crossings
% OUTPUTS
%   -rippleStats [rejected ripples removed]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general parameters
ripple_band = [70 100];

% spike detection parameters
patch_param.hf_shift = 0.04*sfreq;
patch_param.spike_pad = 0*sfreq;
patch_param.spike_score_cutoff = 130;
patch_param.spike_factor = 5/(sfreq/500);
patch_param.hf_factor = 13;

% offset so that indices from detectSpikes align with sharp component
spike_offset = 30 * (sfreq/1000);

% reject if ratio of the largest to third largest peak-to-peak exceeds 3
patch_param.peak_ratio_thresh = 2.5;
if any(isnan(LFPdat))
    nanMask = isnan(LFPdat); 
    LFPdat(nanMask) = 0; 
else
    nanMask = false(1,length(LFPdat));
end

% bandpass at ripple band
[b,a] = butter(3,ripple_band/(sfreq/2));
dat.ripple_band = NaN(size(LFPdat));
for ch = 1:size(dat.ripple_band,1)
    dat.ripple_band(ch,:) = filtfilt(b,a,LFPdat(ch,:));
end

if do_plots
    figure;
end

%%

% NC + HC data
LFPdat_all = [LFPdat; LFPdat_other];
if checkSpike
% check if any overlapping IIS
    spike_mask = detectSpikes_IAV(LFPdat_all', sfreq, [], patch_param.hf_factor, patch_param.spike_factor,...
        patch_param.spike_score_cutoff, patch_param.spike_pad, patch_param.hf_shift);
else
    spike_mask = false(1,length(LFPdat_all));
end
LFPdat(nanMask) = nan;
spike_mask = spike_mask';
spike_inds = find(spike_mask) - spike_offset;

% spikes on any channel
spike_inds_all = find(any(spike_mask,1));

% compute initial average density
pre_density = mean(cellfun(@mean, rippleStats.density));

for ch = 1:numel(rippleStats.locs)
    
    spike_inds = find(spike_mask(ch,:)) - spike_offset;
    
    % 1 to reject
    rej_ind = false(numel(rippleStats.locs{ch}),1);
    
    for rip = 1:numel(rippleStats.locs{ch})
        rej_reason = 'none';
        ampScore = NaN;
        % extract ripple data
        rip_window = rippleStats.window{ch}(rip,1):rippleStats.window{ch}(rip,2);
        
        % move to else clause once done debugging
        ripBand_dat = dat.ripple_band(ch,rip_window);
        ripLFP_dat = LFPdat(ch,rip_window);
        % find zero crossings of ripple bandpass
        zc = diff(ripBand_dat>0);
        zc_type = zc(zc~=0); % rising or falling phase
        zc_t = find(zc==1 | zc==-1); % zero crossing times
            
        % reject if ripple is within spike_proximity_thresh from an IIS on the same channel
        if sum(ismember(rip_window(1)-round(sfreq/2):rip_window(end)+round(sfreq/2), spike_inds))
            rej_ind(rip) = 1;
            rej_reason = 'same ch spike';
        % reject if ripple coincides with a spike on any channel
        elseif sum(ismember(rip_window, spike_inds_all))
            rej_ind(rip) = 1;
            rej_reason = 'other ch spike';
        % reject if a ripple has a dominant peak-to-peak component (>3:1 first to third largest peak-to-peak)
        else

            % compute amplitudes for each rising and falling phase of the LFP
            amps = NaN(numel(zc_type)-2,1);
            for t = 1:numel(zc_type)-2
                if zc_type(t) == -1
                    amps(t) = max(ripLFP_dat(zc_t(t+1):zc_t(t+2))) - min(ripLFP_dat(zc_t(t):zc_t(t+1)));
                else
                    amps(t) = max(ripLFP_dat(zc_t(t):zc_t(t+1))) - min(ripLFP_dat(zc_t(t+1):zc_t(t+2)));
                end
            end
            ampSorted = sort(amps,'descend');
            
            % ripple has fewer than 3 full peak-to-peaks
            if numel(ampSorted) < 3
                rej_ind(rip) = 1;
                rej_reason = 'too short';
            else
                ampScore = ampSorted(1)/ampSorted(3);

                % plot ripples to be rejected

                if ampScore > patch_param.peak_ratio_thresh
                    rej_ind(rip) = 1;
                    rej_reason = 'big cycle';
                end                
            end
        end
        
        % plot examples with rejection reason if applicable
        if do_plots
            clf
            
            yyaxis left
            plot(ripLFP_dat)
            hold on
            vline(zc_t(zc_type==1), 'k')
            vline(zc_t(zc_type==-1), 'r')
            title([rej_reason ' ' num2str(ampScore)])
            
            yyaxis right
            plot(ripBand_dat); hold on;
            waitforbuttonpress
        end
    end
    
    % note: not modifying centeredInd and locs_sleep
    
    if ~isempty(rej_ind)
        % remove bad indices
        rippleStats.HGz{ch}(rej_ind) = [];
        rippleStats.oscFreq{ch}(rej_ind) = [];
        rippleStats.rippleAmp{ch}(rej_ind) = [];
        rippleStats.window{ch}(rej_ind,:) = [];
        rippleStats.locs{ch}(rej_ind) = [];
        rippleStats.duration{ch}(rej_ind) = [];

        % re-calculate
        rippleStats.InterRipPeriod{ch} = diff(rippleStats.locs{ch});
    
    end
    
    rippleStats.density{ch} = sfreq * 60 * numel(rippleStats.locs{ch}) / length(LFPdat);
    
    % could be later cross-referenced to backed up rippleStats if needed
    rippleStats.patch{ch} = rej_ind;
    
end

% add patch metadata
rippleStats.patch_param = patch_param;
rippleStats.patch_param.density_change = mean(cellfun(@mean, rippleStats.density)) - pre_density;

% save rippleStats in a new folder here or return it first

end