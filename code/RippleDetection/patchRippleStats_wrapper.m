function patchRippleStats_wrapper(subjName)
%%



% if using as a script
%subject_list = {'CC08'}%subjects(11:18);%{'CC60'}%'CC30','CC31','CC39','CC49','CC55','CC60','CC69','CC91'};
% subject_list = {'MG67_UA_notch'};

if contains(subjName,'CC')
    % CC patients have SEEG NC and HC

    subjects = {'CC04','CC08','CC15','CC17','CC18','CC20','CC23','CC24','CC25','CC26',...
    'CC31','CC39','CC49','CC55','CC60','CC69','CC91'};

    [~,subjs] = intersect(subjects, subject_list);

    % all sleep files
    sleep_files = {[2,3,7],[2,5,9,10],[3,5,6,7],[1],...
            [2,3,5],[4,10,12],[5,6,7,9,11],[2,3,6,8],[1,2,4,5,7],...
            [3,7,8],[1,2,5],[2,4,5,6,7,8,10,11],[3,8,9,12],...
            [2,6,8],[2,6,7],[1,2,3,5,7,10],[1,3,4,5]};
else
    subjects = {subjName};
    subjs = 1;
    sleep_files = {1};
end

% subject = 'CC08';
conds = {'waking'};%{'waking', 'sleep'};

sfreq = 1000;

do_plots = 1;

for s = 1:numel(subjs)
    subj = subjs(s);
    subject = subjects{subj}
    for c = 1:numel(conds)
        cond = conds{c}
        tic
        
        clear data sleepdata dat rippleStats rippleStats_HC rippleStats_NC
    
        % load labels
        if contains(subject,'CC')
            load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/labels/%s_labels.mat', subject), 'labels');
        end
        % load NC and HC data
        dat.raw_NC = [];
        dat.raw_HC = [];
        if strcmp(cond, 'sleep')
            longdata_dir = sprintf('/space/seh8/1/halgdev/projects/xjiang3/Ripples_newSubj/%s/longdata',subject);    
            load(sprintf('%s/longdata_%s_%d.mat',longdata_dir,subject,sleep_files{subj}(1)),'ctx_ind')
            for sleep_ind = 1:numel(sleep_files{subj})
                sleep_ID = sleep_files{subj}(sleep_ind);
                load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/segment_indices/sleep_%d_N23_seg_NC.mat', subject, sleep_ID))
                dat.raw_NC = [dat.raw_NC sleepdata]; clear sleepdata
                load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/noWaveletReject/HC/sleep_sleep%d_N23.mat', subject, sleep_ID));
                dat.raw_HC = [dat.raw_HC sleepdata]; clear sleepdata
            end

        elseif strcmp(cond, 'wake')
            if contains(subject,'CC')
                load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/waking/HC/waking_1.mat', subject))
                dat.raw_HC = data; clear data
                load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/waking/NC/waking_1.mat', subject))
                dat.raw_NC = data; clear data
                clear data
            else
                load(sprintf('/space/seh9/5/halgdev/projects/iverzh/ripple/NC_ripple/%s/waking/ALLCHAN/waking_1.mat', subject))
                dat.raw_NC = data; clear data
            end
        end
        
        %% NC
            
        % load NC rippleStats
        load(sprintf('/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles/%s_ripple_stats_%s_ALLCHAN_.mat', subject, cond), 'rippleStats')
        rippleStats_NC = rippleStats; clear rippleStats;

        if sum(rippleStats_NC.recordingLength) ~= size(dat.raw_NC,2)
            error('NC rippleStats does not match NC data length')
        end

        % run NC chanSelect
        if contains(subject,'CC')
            orig_chan_labels = rippleStats_NC.chanLabels;
            [rippleStats_NC, chan_labels_NC, ~, dat.raw_NC] = chanSelect(subject, rippleStats_NC, labels, dat.raw_NC);
            ephys_ch_NC = ismember(orig_chan_labels, chan_labels_NC); clear orig_chan_labels
        end
        % patch and save NC rippleStats
        rippleStats = patchRippleStats(rippleStats_NC, dat.raw_NC, dat.raw_HC, sfreq, do_plots);
        if contains(subject,'CC')
            rippleStats.ephys_ch = ephys_ch_NC; % mark which ephys channels were removed
        end
        save(sprintf('/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles/%s_ripple_stats_%s_patched.mat', subject, cond), 'rippleStats');
        clear rippleStats

        %% HC

        % load HC rippleStats
        load(sprintf('/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles/%s_ripple_stats_%s_HC_highPrec.mat', subject, cond), 'rippleStats')
        rippleStats_HC = rippleStats; clear rippleStats;

        if sum(rippleStats_HC.recordingLength) ~= size(dat.raw_HC,2)
            error('HC rippleStats does not match NC data length')
        end

        % run HC chanSelect
        orig_chan_labels = rippleStats_HC.chanLabels;
        [rippleStats_HC, chan_labels_HC, ~, dat.raw_HC] = chanSelect(subject, rippleStats_HC, labels, dat.raw_HC);
        ephys_ch_HC = ismember(orig_chan_labels, chan_labels_HC); clear orig_chan_labels
        
        % patch and save HC rippleStats
        rippleStats = patchRippleStats(rippleStats_HC, dat.raw_HC, dat.raw_NC, sfreq, do_plots);
        rippleStats.ephys_ch = ephys_ch_HC; % mark which ephys channels were removed
        save(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/rippleStats/%s_ripple_stats_%s_HC_highPrecision.mat', subject, cond), 'rippleStats');
        clear rippleStats
%         
        toc
    end
end




