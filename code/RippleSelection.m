% RippleSelection.m 
% Script used to detect ripple oscillations in SEEG data. Methods are described
% in the following manuscript:
% 
% Widespread ripples synchronize human cortical activity during sleep, waking, 
% and memory recall. CW Dickey, IA Verzhbinsky, X Jiang, BQ Rosen, S Kajfez, 
% B Stedelin, JJ Shih, S Ben-Haim, AM Raslan, EN Eskander, J Gonzalez-Martinez,
% SS Cash, E Halgren
% 
% Ilya A. Verzhbinsky, Halgren Lab, 02.24.2022
clc
clear 
close all

addpath(genpath('.'))


%% Inputs



% sEEG subjects 
subj_list_full = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
                  'S11','S12','S13','S14','S15', 'S16','S17'};
segmentFiles_set = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; %list of segement Ids used for each segment (if multiple segments used).

runList = 1;%:length(subj_list_full); 
recordingState = 'sleep'; % sleep or wake
modifier = '';
IISflag = 1; %1 - check for spikes. 0 - do not check for spikes

detectRipples = true; %true - perform preliminary ripple detection. false - load previous prelim detection
saveRippleStats = true; %save selected ripples

fs = 1000; %sampling rate
rippleband = [70 100]; %ripple band in Hz.

rippPrelimDir = './'; % If detectRipples is false, path to preliminary ripple detection mat file. 
                      % If detectRipples is true, path to preliminary ripple detection output folder.
LFPdir ='./'; % path to LFP mat file. LFP variable should be named 'data' and be formatted as nChannels x time.
if strcmp(recordingState,'sleep')
    isSleep = 1;
elseif strcmp(recordingState,'waking')
    isSleep = 0;
else
    error("incorrect RecordingState argument. Can only be 'sleep', 'waking'")
end

masterFolder = './../data/Figures';
if ~isfolder(masterFolder); mkdir(masterFolder); end

detectExportFolder = './../data/prelimDetection';
if ~isfolder(detectExportFolder); mkdir(detectExportFolder); end

matExportFolder = './../data/matFiles';
if ~isfolder(matExportFolder); mkdir(matExportFolder); end

preProcExportFolder = './../data/preProc';
if ~isfolder(preProcExportFolder); mkdir(preProcExportFolder); end

winSec = 2; % pre and post event window in seconds
spikeCheck = 500; %+/- for spike check (10 Hz) abs(LFP) > 1000
rippleWindowSec = 0.100; %+/- window used to compute ripple ocsillation freq. in seconds
IIScheckSec = 0.500; %+/- in secs

% Rejection Parameters.
% 1 in position 2 indiciates rejejction criterion is active.
% 0 in position 2 indiciates rejejction criterion is not used.

RejectParams.RBHFzscore = [0, 0]; %Reject zscore(Ripple) - zscore(HF) < RBHFzscore
RejectParams.sharpzscore = [7, 1]; %100 Hz highpass filter z score threshold
RejectParams.LFPdiff = [50, 0]; %LFP derivative absolute threshold in microV
RejectParams.LFPdiffzscore = [4, 0]; %LFP derivative z score threshold
RejectParams.RBzscore = [3, 1]; % rippleband z score threshold
RejectParams.LFPthresh = [1000, 1]; %LFP rejection threshold
RejectParams.minDuration = [0.025, 1]; % minimum duration in seconds





%%
fprintf('Running Ripple Selection ...\n')

for subj = runList
    subject = subj_list_full{subj};
    
    fprintf('Loading Subject %s ...\n', subject)

    tag = [recordingState,'_',modifier];
    rippleStatsAllSleepSet = struct([]);
    shift = 0;
    for s = 1:length(segmentFiles_set{subj})
        segment_ID = segmentFiles_set{subj}(s);
            
        if exist('data','var')
            clear data
        end
         
        if exist('sleepdata','var')
            clear sleepdata
        end

        fprintf('Loading sleep segment %i ...\n', segment_ID)
        load(LFPdir) %this may need to be editted to refelct user file architecture
        
        if detectRipples
            load('./RippleDetection/detectionMethods_default.mat')
            methods.separation_flag = 0;
            methods.adjacent_merge = 1;
            methods.adjacency = 25; % in sampling points
            % broader range for freq center investigation
            methods.hp_freq_ST = 60;  % in Hz
            methods.lp_freq_ST = 120;
            methods.looseEval = 1;
            methods.gap_thresh = 0.1;   % set as 100 ms for now; might want to turn this off for NC
            methods.do_pmax = 0;   % flag off since no hand-marked template for cortical ripples (there's no "sharp-wave" equivalent as far as we know)
            methods.do_wavelet = 0;  % next 7 params are to tweak outlier_wavedecompC_dbtest2.m
            methods.time_bin = 2;
            methods.wavthresh_Haar = 6;
            methods.wavthresh_RL1 = 6;
            methods.wavthresh_RL2 = 8;
            methods.wavthresh_RL3 = 8;
            methods.wavthresh_RL4 = 5;
            methods.wavthresh_RL5 = 4;
            methods.HP_reject = 1;
            methods.HP_thresh = 3;   % reject based on high gamma powe
            methods.parpool_size = 16;
% 
            gMsk = ones(1, size(LFP,1));    
            mask = ones(1, size(LFP,2));    
            rippleDetect(data,fs,savedir,1,gMsk, methods, mask, recordingState);
        else

            try
                load(detectExportFolder)
            catch
                WarningString = ['could not load ',f];
                warning(WarningString)
                saveRippleStats = false;
                break
            end
        end %detect prelim ripples

        if exist('sfreq','var')
            fs = sfreq;
        end
        
        if exist('all_ripple_locs','var')
            good_ripple_locs_set = all_ripple_locs;
        end

        if exist('data','var')
            BroadbandData = data;
            clear data
        elseif exist('sleepdata','var')
            BroadbandData = sleepdata;
            clear sleepdata
    
        end
        
       
        try 
            load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/%s/nan_edge_mask.mat', subject, recordingState),'nan_edge_mask') 
        catch 
            nan_edge_mask = false(1,size(BroadbandData,2));
        end
        

        if strcmp(recordingState,'sleep')
            polCheckFiles = dir(fullfile(preProcExportFolder,[subject,'_',num2str(segment_ID),'_',tag,'_polCheck.mat']));

            if isempty(polCheckFiles)
                polCheck = checkPolarity(BroadbandData,fs); %check polarity by comparing slow oscillations w HG
                try
                    save(fullfile(preProcExportFolder,[subject,'_',num2str(segment_ID),'_',tag,'_polCheck.mat']), 'polCheck', '-v7.3')
                catch
                    warning('could not save polarity check')
                end
            else
                load(fullfile(preProcExportFolder,[subject,'_',num2str(segment_ID),'_',tag,'_polCheck.mat']))
            end
        else
            polCheck = ones(size(BroadbandData));
        end

        BroadbandData = polCheck .* BroadbandData;

        if IISflag
            spikeCheckFiles = dir(fullfile(preProcExportFolder,[subject,'_',num2str(segment_ID),'_',tag,'_spikeCheck.mat']));
            if isempty(spikeCheckFiles)
                [spikeMask,spikeBounds] = detectSpikes_IAV(BroadbandData', fs);
                try
                    save(fullfile(preProcExportFolder,[subject,'_',num2str(segment_ID),'_',tag,'_spikeCheck.mat']),'spikeBounds','spikeMask', '-v7.3')
                catch
                    warning('could not save spike check')
                end
            else
                load(fullfile(preProcExportFolder,[subject,'_',num2str(segment_ID),'_',tag,'_spikeCheck.mat']))
            end
        else 
            spikeMask = zeros(size(BroadbandData));
        end

        if isa(chan_labels,'double'); chan_labels = cellstr(string(chan_labels)); end

        for ch = 1:3 %size(BroadbandData,1)
            initialRipples = good_ripple_locs_set{ch};
   
            fsNew = fs;
            win = round(winSec * fsNew); %win in samples.
            rippleWindow = round(rippleWindowSec * fsNew); %ripple window in samples
            IIScheck = round(IIScheckSec * fsNew); %spike check window in samples 
            
            checkDiff = [diff(initialRipples) 5];
            initialRipples(checkDiff<5) = [];
            [rippleSleepSet,rejectVec] = AnalyzeRipple(BroadbandData, rippleband, win, rippleWindow, ...
                                                       IIScheck, spikeMask, nan_edge_mask, IISflag, ...
                                                       ch, initialRipples, RejectParams, fsNew, ...
                                                       chan_labels, shift,1);

            rippleStatsSleepSet =  computeRippleStats(rippleSleepSet, fs);




            rippleStatsSleepSet.locs{1} = rippleSleepSet.goodRipples;
            rippleStatsSleepSet.rejectVec{1} = rejectVec;
            density = length(rippleStatsSleepSet.locs{1} / length(BroadbandData)) * fs * 60;
            rippleStatsSleepSet.density{ch} = density;
            
            rippleStatsSleepSet = patchRippleStats(rippleStatsSleepSet, BroadbandData(ch,:), [], fs, 0);

            if isempty(rippleStatsAllSleepSet)
                fn = fieldnames(rippleStatsSleepSet);
                for ch2 = 1:length(chan_labels)
                    for f = 1:numel(fn)
                        rippleStatsAllSleepSet(1).(fn{f}){ch2} = [];
                    end

                end
            end


            fn = fieldnames(rippleStatsSleepSet);
            if numel(rippleStatsSleepSet.locs{1}) > 5
                for f = 1:numel(fn)
                    if isa(rippleStatsSleepSet.(fn{f}), 'cell')
                        fDat = rippleStatsSleepSet.(fn{f}){1};
                        DIM = size(fDat,1);
                        if DIM > 1 || size(fDat,2) == 4001
                            rippleStatsAllSleepSet.(fn{f}){ch} = [rippleStatsAllSleepSet.(fn{f}){ch}; fDat];
                        elseif DIM == 1
                            rippleStatsAllSleepSet.(fn{f}){ch} = [rippleStatsAllSleepSet.(fn{f}){ch}, fDat];
                        end
                    elseif isa(rippleStatsSleepSet.(fn{f}), 'struct')
                        fDat = rippleStatsSleepSet.(fn{f});
                       rippleStatsAllSleepSet.(fn{f}) = fDat;
                    end

                end
            end

            rippleStatsAllSleepSet.rejectVec{ch} = rejectVec;


        end

        rippleStatsAllSleepSet.recordingLength(s) = size(BroadbandData,2);

        locs_sleep = nan(1, 5e5); %matrix to store ripple locs by sleep set
        for ch = 1:size(BroadbandData,1)
                   %save ripple indicies by sleep set
                locs_sleep(1:length(rippleStatsAllSleepSet.goodRipples{ch})) = rippleStatsAllSleepSet.locs{ch};
                rippleStatsAllSleepSet.locs_sleep{ch}(s,1:size(locs_sleep,2)) = locs_sleep;
        end

        shift = shift + size(BroadbandData,2); %shift for ripple indexing


    end

    %% Exporting Ripple Data
    if saveRippleStats
        plotStats =  computeRippleStats(rippleAllSleepSet, fs);

        rippleStats = plotStats;
        for ch = 1:size(BroadbandData,1)
            rippleStats.locs{ch} = rippleAllSleepSet.goodRipplesConcat{ch};
            rippleStats.locs_sleep{ch} = rippleAllSleepSet.locs_sleep{ch};


            density = length(rippleAllSleepSet.goodRipples{ch}) / sum(rippleAllSleepSet.recordingLength) * fs * 60;
            rippleStats.density{ch} = density;
        end

        rippleStats.recordingLength = rippleAllSleepSet.recordingLength;
        rippleStats.fs = fs;
        rippleStats.chanLabels = chan_labels;
        rippleStats.RejectParams = RejectParams;
        rippleStats.rejectVec = rippleAllSleepSet.rejectVec;
        rippleStats.RB = rippleband;
        rippleStats.IISflag = IISflag;



        save(sprintf('%s/%s_ripple_stats_%s.mat', matExportFolder, subject, tag), 'rippleStats', '-v7.3')
    end

    
    
end
    
    
    
   