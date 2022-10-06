


parpool(20)

%%
clc
clear 
close all

addpath(genpath('/home/iverzh/ripple/CortRipple'))

addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple'))


%% Inputs
subj_list_full = {'CC04'};
%               

sleepfiles_set = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

win=3000;
binSize = 1; %ms
hist_bins = 2*win/binSize;
nIter = 200;

recordingState = 'sleep';

rippleDir = '/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles/';
badChanDir = '/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/bad_chan';

saveDir = ['/space/seh9/2/halgdev/projects/iverzh/ripples/PRTH/',recordingState,'/NCNC_patch_highPrec'];
if ~isfolder(saveDir); mkdir(saveDir); end





%%

c = 0;
for subj = 1:length(subj_list_full)%[1,2,3,4,6,9,10,12,13,16]
    
    subject = subj_list_full{subj};
    
    fprintf([subject, '\n'])
    
    subjPRTH = [];
    
    
    load(fullfile(parcelDir,[subject,'_labels.mat']));
    switch recordingState
        case 'sleep'
            load(fullfile(rippleDir,[subject,'_ripple_stats_sleep.mat']))
        case 'wake'
            load(fullfile(rippleDir,[subject,'_ripple_stats_wake.mat']))
    end
    
    [rippleStats, chanLabels, labelsNC] = chanSelect(subject,rippleStats,labels,[]);


    for ch = 1:length(rippleStats.locs)
        


        for ch2 = 1:length(rippleStats.locs)
            
         
            
            if ch ~= ch2 
                nRipples        = length(rippleStats.locs{ii});
    
                % Compute PRTH
                [event_PRTH, null_PRTH] = comptutePRTH(rippleStats, recordingState, ch, win, binSize, ...
                                            nIter, '', 'NC', sleepfiles_set{subj}, subject, ch2,0);

    
                subjPRTH.eventPRTH{ch,ch2} = event_PRTH;
                subjPRTH.nullPRTH{ch,ch2} = null_PRTH;
                subjPRTH.nRipples{ch,ch2} = nRipples;
    
    
                fprintf([subject, ' ' ,chan,' ',chan2,' ' ])
                fprintf(' done.\n')
                
                c = c + 1;
            end
        
        end

       
    end
    subjPRTH.chanLabels = labelsNC(:,1);
    save(fullfile(saveDir,[subject, '_NCNC_PRTH.mat']), 'subjPRTH');
    fprintf('\n\n')
end



