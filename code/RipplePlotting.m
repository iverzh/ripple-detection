% RipplePlotting.m 
% Script used to plot ripple oscillations detected in SEEG data. Methods are described
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

addpath(genpath('./'))

addpath(genpath('./eeglab14_1_2b/')) %need to download eeglab
rmpath(genpath('./eeglab14_1_2b/functions/octavefunc'))

%% Paths / Inputs


subj_list_full = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10', ...
                  'S11','S12','S13','S14','S15', 'S16', 'S17'};
segmentFiles_set = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; %list of segement Ids used for each segment (if multiple segments used).


% sleepfiles_set = {1,1,1,1,1,1,1,1,1,1,1,1};
runList = 1:length(subj_list_full); 
winSec = 2; % pre and post event window in seconds
rippleWindowSec = 0.100; %+/- window used to compute ripple ocsillation freq. in seconds


recordingState = 'sleep'; %sleep or waking
tag = recordingState;
modifier = '';

exportFolder = './../Figures/';
if ~isfolder(exportFolder); mkdir(exportFolder); end              
              
rippleStatsFolder = './../data/matFiles';              
preProcExportFolder = './../data/preProc';
LFPdir ='./'; % path to LFP mat file. LFP variable should be named 'data' and be formatted as nChannels x time.


%% Load rippleStats and compute bandpasses
for subj = runList %1:length(subjID)
    
    subject = subj_list_full{subj};

    load(fullfile(rippleStatsFolder,[subject,'_ripple_stats_',recordingState,'_',location,'_',modifier,'.mat']));

    data = [];

    segmentFiles = segmentFiles_set{subj};

    for si = 1:length(segmentFiles)
        s = segmentFiles(si);
        load(LFPdir) %this may need to be editted to refelct user file architecture

        if isempty(data)
            data = temp.sleepdata .* polCheck;
        else
            data = [data, temp.sleepdata .* polCheck];
        end
    end



    win = winSec * rippleStats.fs;
    rippleWindow = rippleWindowSec * rippleStats.fs;

    chan_labels = rippleStats.chanLabels;
    rippleband = rippleStats.RB;
    for ch = 1:size(data,1)

        [rippleSleepSet,~] = AnalyzeRipple(data, rippleband, win, rippleWindow, 0, 0, nan_edge_mask, 0, ...
                                                       ch, rippleStats.locs{ch}, [], rippleStats.fs, chan_labels, 0,0);

        LFPbandpasses = [];                                           
        fn = fieldnames(rippleSleepSet);
        %      rippleAllSleepSet = rippleSleepSet;
        for f = 1:numel(fn)
            LFPbandpasses(1).(fn{f}){ch} = [];
        end


        fn = fieldnames(rippleSleepSet);
        if numel(rippleSleepSet.rippleAmp) > 5
            for f = 1:numel(fn)
                fDat = rippleSleepSet.(fn{f});
                DIM = size(fDat,1);
                if DIM > 1 || size(fDat,2) == 4001
                    LFPbandpasses.(fn{f}){ch} = [LFPbandpasses.(fn{f}){ch}; fDat];
                elseif DIM == 1
                    LFPbandpasses.(fn{f}){ch} = [LFPbandpasses.(fn{f}){ch}, fDat];
                end

            end
        end

        LFPbandpasses.recordingLength = size(data,2);

        %% Plotting
        fprintf(sprintf('exporting plots for channel %s\n', chan_labels{ch}))

        summaryFolder = fullfile(exportFolder, subject, ['ChannelSummary_', tag, '_',modifier]);
        singleTrialFolder = fullfile(exportFolder, subject, ['SingleTrials_',tag,'_', modifier], chan_labels{ch});
        TFexportFolder = fullfile(exportFolder, subject, ['TF_', tag,'_',modifier]);


        if ~isfolder(fullfile(summaryFolder,'PNG')); mkdir(fullfile(summaryFolder, 'PNG')); end
        if ~isfolder(fullfile(summaryFolder,'PDF')); mkdir(fullfile(summaryFolder, 'PDF')); end
        if ~isfolder(fullfile(summaryFolder,'Fig')); mkdir(fullfile(summaryFolder, 'Fig')); end

        if ~isfolder(fullfile(singleTrialFolder,'PNG')); mkdir(fullfile(singleTrialFolder, 'PNG')); end
        if ~isfolder(fullfile(singleTrialFolder,'PDF')); mkdir(fullfile(singleTrialFolder, 'PDF')); end
        if ~isfolder(fullfile(singleTrialFolder,'Fig')); mkdir(fullfile(singleTrialFolder, 'Fig')); end

        if ~isfolder(fullfile(TFexportFolder,'PNG')); mkdir(fullfile(TFexportFolder, 'PNG')); end
        if ~isfolder(fullfile(TFexportFolder,'PDF')); mkdir(fullfile(TFexportFolder, 'PDF')); end
        if ~isfolder(fullfile(TFexportFolder,'Fig')); mkdir(fullfile(TFexportFolder, 'Fig')); end




        if numel(LFPbandpasses.rippleAmp{ch}) >= 20
           meanLFP = mean(LFPbandpasses.raw{ch});

            for rip_num = 1:size(LFPbandpasses.raw{ch},1)
               rippleCOV = cov(LFPbandpasses.raw{ch}(rip_num,:), meanLFP);
               LFPbandpasses.cov{ch}(rip_num) = rippleCOV(1,2);     
            end

    %             h = chanMean_SingleTrials(LFPbandpasses, ch, win, fs, subject, chan_labels);
    %             saveas(h, sprintf('%s/%s_%s_summarySingleTrials.png', summaryFolder, subject, chan_labels{ch}))

            hsummary = chanOverviewPlot(LFPbandpasses, rippleStats, ch, win, rippleStats.fs, subject, chan_labels);
            savepdf(hsummary, sprintf('%s/PDF/%s_%s_summary.pdf', summaryFolder, subject, chan_labels{ch}))
            saveas(hsummary,  sprintf('%s/PNG/%s_%s_summary.png', summaryFolder, subject, chan_labels{ch}))
            saveas(hsummary,  sprintf('%s/Fig/%s_%s_summary.fig', summaryFolder, subject, chan_labels{ch}))

    %                     
            hTF = chanTFplot(LFPbandpasses, ch, win, rippleStats.fs, subject, chan_labels,'TF');
            savepdf(hTF, sprintf('%s/PDF/%s_%s_TimeFreq.pdf', TFexportFolder, subject, chan_labels{ch}))
            saveas(hTF,  sprintf('%s/PNG/%s_%s_TimeFreq.png', TFexportFolder, subject, chan_labels{ch}))
            saveas(hTF,  sprintf('%s/Fig/%s_%s_TimeFreq.fig', TFexportFolder, subject, chan_labels{ch}))
    % % % 
            hTF = chanTFplot(LFPbandpasses, ch, win, rippleStats.fs, subject, chan_labels,'LFP');
            savepdf(hTF, sprintf('%s/PDF/%s_%s_MeanLFP.pdf', TFexportFolder, subject, chan_labels{ch}))
            saveas(hTF,  sprintf('%s/PNG/%s_%s_MeanLFP.png', TFexportFolder, subject, chan_labels{ch}))
            saveas(hTF,  sprintf('%s/Fig/%s_%s_MeanLFP.fig', TFexportFolder, subject, chan_labels{ch}))

          %  Plot by rippleBand amp 
            [~,I] = sort(LFPbandpasses.rippleAmp{ch}, 'descend');
            for f = 1:5
                h3 = figure('Position', [3 30 1908 885], 'Visible', 'off');
                for i = 1:4

                    rip_num = I(4*(f-1) + i);
                    rippleFullPlot(LFPbandpasses, i, ch, rip_num, rippleStats.fs, rippleband)

                end
                orient(h3,'landscape')
                savepdf(h3, sprintf('%s/PDF/%s_%s_exampleRipple_%i.pdf', singleTrialFolder, subject, chan_labels{ch}, f))
                saveas(h3,  sprintf('%s/PNG/%s_%s_exampleRipple_%i.png', singleTrialFolder, subject, chan_labels{ch}, f))
                saveas(h3,  sprintf('%s/Fig/%s_%s_exampleRipple_%i.fig', singleTrialFolder, subject, chan_labels{ch}, f))
                close
            end


            close all
        end
    end
end