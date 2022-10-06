

% 
% INPUTS:
%     - rippleStats (NC rippleStats structure)
%     - recordingState - 'sleep' or 'wake';
%     - ch 
%     - win: +- window in ms (default 3000)
%     - binWidth: in ms (default 1)
%     - nIter
%     - mode: 
%     - sleepfiles_set: needs to be single integer if 'wake'
%     - flipFlag: 0 - NC at t = 0, 1 - flipped.



function [event_PRTH, null_PRTH] = computePRTH(rippleStats, recordingState, ch, win, binWidth, nIter, ...
                                                histEventsPath, mode, sleepfiles_set, subject, ch2, ...
                                                flipFlag)


edges = -win+(binWidth/2):binWidth:win;

peri_event = nan(1e8,1);
null_PRTH  = zeros(nIter,length(edges)-1);

count = 1;
for s = 1:length(sleepfiles_set)
    sleep_ID = sleepfiles_set(s);
    
    switch recordingState
        case 'sleep'
            try
                load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/N23_masks/%s_%02d_N23_mask.mat', subject, sleep_ID));     
            catch
                warning('off','backtrace')
                warning('Cannot load NREM mask')
                NREM_mask = ones(1,round(1000 / rippleStats.fs * rippleStats.recordingLength(s)));
            end
        case 'wake'
            NREM_mask = ones(1,round(1000 / rippleStats.fs * rippleStats.recordingLength(s)));
    end
                
    NREM_mask_ind = find(NREM_mask);
    
    % remove ripples with indices less than window size
    ripple_inds = rippleStats.locs{ch};
    ripple_inds(isnan(ripple_inds)) = [];
 
    if s == 1
        ripple_inds = ripple_inds(ripple_inds <= sum(rippleStats.recordingLength(1))); 
    else
        keep = ripple_inds <= sum(rippleStats.recordingLength(1:s)) & ripple_inds > sum(rippleStats.recordingLength(1:s-1));
        ripple_inds = ripple_inds(keep);
        ripple_inds = ripple_inds - sum(rippleStats.recordingLength(1:s-1));
    end
    
    ripple_inds(ripple_inds<win) = [];


    
    switch mode
        case 'US'
            if contains(histEventsPath,'xjiang')
                SOFilename = sprintf('%s/%s_sleep%i_SO.mat', ...
                                        histEventsPath, subject, sleep_ID);
                load(SOFilename)
                goodCh = find(rippleStats.ephys_ch);
                chEvent = goodCh(ch);
            else
                SOFilename = sprintf('%s/%s_N23_SO.mat', ...
                                        histEventsPath, subject);
                load(SOFilename)
                chEvent = find(contains(chlab,rippleStats.chanLabels{ch}));
                
            end
                
            if ~isempty(chEvent)
                % select channel and down-states (negative peaks only)
                hist_events = SOmat(SOmat(:,3)==chEvent & SOmat(:,2) > 0,1);
            else
                hist_events = [];
            end
            
        case 'DS'
            if contains(histEventsPath,'xjiang')
                SOFilename = sprintf('%s/%s_sleep%i_SO.mat', ...
                                        histEventsPath, subject, sleep_ID);
                load(SOFilename)
                goodCh = find(rippleStats.ephys_ch);
                chEvent = goodCh(ch);
            else
                SOFilename = sprintf('%s/%s_N23_SO.mat', ...
                                        histEventsPath, subject);
                load(SOFilename)
                chEvent = find(contains(chlab,rippleStats.chanLabels{ch}));

            end
            
            if ~isempty(chEvent)
                % select channel and down-states (negative peaks only)
                hist_events = SOmat(SOmat(:,3)==chEvent & SOmat(:,2) < 0,1);
            else
                hist_events = [];
            end
            
       case 'spin'
           if contains(histEventsPath,'xjiang')
                spinFilename = sprintf('%s/spindle_results_%s_Sleep%i_SPINmat.mat', ...
                                        histEventsPath, subject, sleep_ID);
                load(spinFilename)
                goodCh = find(rippleStats.ephys_ch);
                chEvent = goodCh(ch);
           else
               spinFilename = sprintf('%s/%s_SPINmat.mat', ...
                                        histEventsPath, subject);
               load(spinFilename)

               chEvent = find(contains(chlab,rippleStats.chanLabels{ch}));

           end
           
           if ~isempty(chEvent)
               hist_events = spin_locs{chEvent,1};
           else
               hist_events = [];
           end
            
        case 'AMY'
            TBFilename = sprintf('%s/spindle_data_out_%s_sleep%i_TH.mat', ...
                                    histEventsPath, subject, sleep_ID);
            AMY = load(TBFilename);
            
            hist_events = AMY.rippleStats.locs{ch2};
            
        case 'HC'
%             HC = load(fullfile(histEventsPath,[subject,'_ripple_stats_HC_NSWSSR.mat']));
            
            switch recordingState
                case 'sleep'
                    HC = load(fullfile(histEventsPath,[subject,'_ripple_stats_sleep_HC.mat']));
                case 'wake'
                    HC = load(fullfile(histEventsPath,[subject,'_ripple_stats_wake_HC.mat']));
            end
            
            [rippleStatsHC, ~, ~] = chanSelect(subject,HC.rippleStats,'',[]);
            
            hist_events = rippleStatsHC.locs{ch2};

%             count2 = 1;
%             for chHC = 1:length(chanLabelsHC)
%                 iiHC = rippleStatsHC.locs{chHC};
%                 iiHC(isnan(iiHC)) = [];
%                 hist_events(count2 : count2 + length(iiHC) - 1) = iiHC;
%                 count2 = count2 + length(iiHC);
%             end
            
            
%             hist_events = HC.rippleStatsNSWSSR.locs{1};
            hist_events(hist_events == 0) = [];
            
            hist_events = sort(hist_events);

            
            if s == 1
                hist_events = hist_events(hist_events <= sum(HC.rippleStats.recordingLength(1))); 
            else
                keep = hist_events <= sum(HC.rippleStats.recordingLength(1:s)) & hist_events > sum(HC.rippleStats.recordingLength(1:s-1));
                hist_events = hist_events(keep);
                hist_events = hist_events - sum(HC.rippleStats.recordingLength(1:s-1));
            end
            
            hist_events = NREM_mask_ind(hist_events);
            
             if rippleStats.recordingLength ~= rippleStatsHC.recordingLength
                error('recording lengths do not match')
             end
             
        case {'SWR','SSR'}
%             HC = load(fullfile(histEventsPath,[subject,'_ripple_stats_HC_NSWSSR.mat']));
            
            HC = load(fullfile(histEventsPath,[subject,'_ripple_stats_HC_',mode,'.mat']));
            
            structname = ['rippleStats',mode];
            [rippleStatsHC, ~, ~] = chanSelect(subject,HC.(structname),'',[]);
            
            hist_events = rippleStatsHC.locs{ch2};

%             count2 = 1;
%             for chHC = 1:length(chanLabelsHC)
%                 iiHC = rippleStatsHC.locs{chHC};
%                 iiHC(isnan(iiHC)) = [];
%                 hist_events(count2 : count2 + length(iiHC) - 1) = iiHC;
%                 count2 = count2 + length(iiHC);
%             end
            
            
%             hist_events = HC.rippleStatsNSWSSR.locs{1};
            hist_events(hist_events == 0) = [];
            
            hist_events = sort(hist_events);

            
            if s == 1
                hist_events = hist_events(hist_events <= sum(HC.(structname).recordingLength(1))); 
            else
                keep = hist_events <= sum(HC.(structname).recordingLength(1:s)) & hist_events > sum(HC.(structname).recordingLength(1:s-1));
                hist_events = hist_events(keep);
                hist_events = hist_events - sum(HC.(structname).recordingLength(1:s-1));
            end
            
            hist_events = NREM_mask_ind(hist_events);
            
             if rippleStats.recordingLength ~= rippleStatsHC.recordingLength
                error('recording lengths do not match')
            end
        
        case 'TH'
%             HC = load(fullfile(histEventsPath,[subject,'_ripple_stats_HC_NSWSSR.mat']));
            TH = load(fullfile(histEventsPath,[subject,'_ripple_stats_wake_TH.mat']));
            [rippleStatsTH, ~, ~] = chanSelect(subject,TH.rippleStats,'',[]);
            
            if rippleStats.recordingLength ~= rippleStatsTH.recordingLength
                error('recording lengths do not match')
            end

            hist_events = zeros(1e8,1);
            
            count2 = 1;
            for chTH = 1:length(chanLabelsTH)
                iiTH = rippleStatsTH.locs{chTH};
                iiTH(isnan(iiTH)) = [];
                hist_events(count2 : count2 + length(iiTH) - 1) = iiTH;
                count2 = count2 + length(iiTH);
            end
            
            
%             hist_events = HC.rippleStatsNSWSSR.locs{1};
            hist_events(hist_events == 0) = [];
            
            hist_events = sort(hist_events);

            
            if s == 1
                hist_events = hist_events(hist_events <= sum(TH.rippleStats.recordingLength(1))); 
            else
                keep = hist_events <= sum(TH.rippleStats.recordingLength(1:s)) & hist_events > sum(TH.rippleStats.recordingLength(1:s-1));
                hist_events = hist_events(keep);
                hist_events = hist_events - sum(TH.rippleStats.recordingLength(1:s-1));
            end
            
            hist_events = NREM_mask_ind(hist_events);
            
            
          
        case 'NC'
%             NC = load(fullfile(histEventsPath,[subject,'_ripple_stats_sleep.mat']));
            hist_events = rippleStats.locs{ch2};
            
            if s == 1
                hist_events = hist_events(hist_events <= sum(rippleStats.recordingLength(1))); 
            else
                keep = hist_events <= sum(rippleStats.recordingLength(1:s)) & hist_events > sum(rippleStats.recordingLength(1:s-1));
                hist_events = hist_events(keep);
                hist_events = hist_events - sum(rippleStats.recordingLength(1:s-1));
            end
            
            hist_events = NREM_mask_ind(hist_events);
            
%         case 'motif'
%             M = load(histEventsPath);
%             motifs = M.RF_cache_SigWake;
%             hist_events = [];
%             for m = 1:numel(motifs{3})
%                 ch_match = find(motifs{3}{m}(:,2) == ch);
%                 if ~isempty(ch_match)
%                     hist_events = [hist_events; motifs{3}{m}(ch_match,1)];
%                 end
%             end
            
         
    end
    
    if ~isempty(hist_events) && ~isempty(ripple_inds) 
        
        hist_events(hist_events < 0) = [];
        hist_events(isnan(hist_events)) = [];
        
%         hist_events = round(1000 / rippleStats.fs * hist_events);
%         ripple_inds = round(1000 / rippleStats.fs * ripple_inds);

        % create mask for spindle onsets
        event_mask = zeros(1e9,1);
        
        if flipFlag %flip PRTH axes
           [~,temp,~] = intersect(NREM_mask_ind,hist_events);
           hist_events = NREM_mask_ind(ripple_inds);
           ripple_inds = temp;
        end

        % extract event times and load into mask
        for event = 1:numel(hist_events)
            event_mask(hist_events(event)) = 1;
        end

        % truncate mask so that it just extends to the end of the spindles
        % this could be done in a better way using the limits of the actual data
        event_mask(max(NREM_mask_ind(ripple_inds))+win+1:end) = [];
        
        ripple_inds(NREM_mask_ind(ripple_inds) < win) = [];


        nullDistributionMat = computeNullDistributionParPool(event_mask, ripple_inds, nIter, NREM_mask_ind, win, edges, rippleStats.fs);
        null_PRTH = null_PRTH + nullDistributionMat;
        
        % generate structure for peri-ripple time histogram of spindle onsets
        for ripple = 1:numel(ripple_inds)
            PRTtemp = find(event_mask(NREM_mask_ind(ripple_inds(ripple))-win:NREM_mask_ind(ripple_inds(ripple))+win))-win;
            if ~isempty(PRTtemp)
                peri_event(count : count + length(PRTtemp) - 1) = PRTtemp;
                count = count + length(PRTtemp);
            end
        end
        
    end
end

peri_event(isnan(peri_event)) = [];

% if any(peri_event > win

% [event_PRTH,~, ~] = histcounts(peri_event(peri_event < win & peri_event > -win), binWidth);
[event_PRTH,~, ~] = histcounts(peri_event, edges);

% event_PRTH = peri_event;
