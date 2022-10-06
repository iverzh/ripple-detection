

parpool(16)
%%

clear
close all
clc
% 1 - frontal, 2 - temporal, 3 - parietal, 4 - occipital, 5 -
% cingulate/limbic, 6- other



addpath(genpath('/home/iverzh/ripple/CortRipple'))

addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple'))



%% Inputs
subj_list_full = {'CC08'};

           

win=3000;
binSize = 1; %25; %ms
hist_bins = 2*win/binSize;
nIter = 100;

recordingState = 'wake';
mode = 'HC';
flipFlag = 1; % 0 - NC-R is on x axis, 1 - HC-R is on x axis 

sleepfiles_set= {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

rippleDir = '/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles/';
% rippleDir = '/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/rippleStats';
rippleDirHC = '/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles/';

saveDir = fullfile('/space/seh9/2/halgdev/projects/iverzh/ripples/PRTH/',recordingState,[mode,'_patch']);
if ~isfolder(saveDir); mkdir(saveDir); end

% badChanDir = '/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/bad_chan/';


%% Stats per parcel

for subj = 1:length(subj_list_full) %5:length(subj_list_full)%[1,2,3,4,6,9,10,12,13,16]
   
    subject = subj_list_full{subj};
    
    fprintf([subject, '\n'])
    
    subjPRTH = [];
     
    
    switch recordingState
        case 'sleep'
            load(fullfile(rippleDir,[subject,'_ripple_stats_sleep.mat']))
%             HC = load(fullfile(rippleDirHC,[subject,'_ripple_stats_sleep_HC.mat']));
            HC = load(fullfile(rippleDirHC,[subject,'_ripple_stats_HC_',mode,'.mat']));
        case 'wake'
            load(fullfile(rippleDir,[subject,'_ripple_stats_wake.mat']))
            HC = load(fullfile(rippleDir,[subject,'_ripple_stats_wake_HC.mat']));
    end

    for ch = 1:length(rippleStats.locs)
        

        nRipples        = length(rippleStats.locs{ch});
        recordingLength = sum(rippleStats.recordingLength);

        
        for chHC = 1:length(HC.rippleStats.locs)
        % Compute PRTH
            HCdir = fullfile(rippleDirHC);%,'..','matFiles');
            [event_PRTH, null_PRTH] = computePRTH(rippleStats, recordingState, ch, win, binSize, nIter, HCdir, mode, ...
                                                   sleepfiles_set{subj}, subject, chHC,flipFlag);

            
            subjPRTH.(mode).nullPRTH{ch,chHC}  = null_PRTH;
            subjPRTH.(mode).eventPRTH{ch,chHC} = event_PRTH;
            subjPRTH.(mode).nRipples{ch,chHC} = nRipples;

        end
        
%                 checkSig = mask2bounds(h(30-10:30+10));
%                 checkSig = abs(checkSig(:,2) - checkSig(:,1)) + 1;
%                 if max(checkSig) > 2
%                     PRTH.(parcelLabel{1}).(mode).nSigChan = PRTH.(parcelLabel{1}).(mode).nSigChan + 1;
%                 end
    %    


        fprintf([chan,' '])
        fprintf('done.\n')
            
    end
        
    
    subjPRTH.chanLabels = rippleStats.chanLabels;
    subjPRTH.chanLabelsHC = HC.rippleStats.chanLabels;
    if ~flipFlag
        save(fullfile(saveDir,[subject, '_',mode,'_PRTH.mat']), 'subjPRTH');
    elseif flipFlag 
         save(fullfile(saveDir,[subject, '_',mode,'_PRTH_flip.mat']), 'subjPRTH');

    end
    fprintf('\n\n')
end

%  save('PRTH.mat','PRTH','parcelstats')



%%

% 




