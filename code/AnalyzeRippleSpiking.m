% AnalyzeRippleSpiking.m 
% Script used to produce figures 5B and 5C in the following manuscript:
% 
% Widespread ripples synchronize human cortical activity during sleep, waking, 
% and memory recall. CW Dickey, IA Verzhbinsky, X Jiang, BQ Rosen, S Kajfez, 
% B Stedelin, JJ Shih, S Ben-Haim, AM Raslan, EN Eskander, J Gonzalez-Martinez,
% SS Cash, E Halgren
% 
% Charles Dickey, Halgren Lab, 02.24.2022

clc
clear 
close all


%% Inputs

dataDir = './../data/UtahArray';
unit_types={'pyr', 'int'};
savedir ='./../data/Figures';
if ~isfolder(savedir); mkdir(savedir); end

subjects = {'U1', 'U2', 'U3'};



%% stats of spiking during ripples

% create plot of all circular means for PY and IN
% stats for individual PY-IN pairs


min_n = 30; % minimum number of spikes

binEdges = 0:(pi/16):(2*pi);


stats.allSig = [];
stats.num_IN_min30spikes = NaN(numel(subjects),1);
stats.num_PY_min30spikes = NaN(numel(subjects),1);

PY_mean_phases = [];
IN_mean_phases = [];
stats.allSig = [];

p_PY = [];
p_IN = [];

for subj = 1:numel(subjects)
    subject = subjects{subj};
    PY = load(sprintf('%s/%s/ripple_phase_%s_pyr.mat', dataDir, subject, 'DR'));
    IN = load(sprintf('%s/%s/ripple_phase_%s_int.mat', dataDir, subject, 'DR')); % FIX THIS TO LOAD DIFFERENT SUBJECTS
    
    stats.sig{subj} = [];
    IN_min_mask = cellfun('size', IN.ripple_phase, 2) >= 30;
    PY_min_mask = cellfun('size', PY.ripple_phase, 2) >= 30;
    stats.num_IN_min30spikes(subj,1) = sum(IN_min_mask);
    stats.num_PY_min30spikes(subj,1) = sum(PY_min_mask);
    stats.num_IN_min30spikes(subj,2) = numel(IN.ripple_phase); % total
    stats.num_PY_min30spikes(subj,2) = numel(PY.ripple_phase);
    IN.ripple_phase = IN.ripple_phase(IN_min_mask);
    PY.ripple_phase = PY.ripple_phase(PY_min_mask);
    
    % compute circular mean phases
    for un = 1:numel(PY.ripple_phase)
        if ~isempty(PY.ripple_phase{un})
            PY_mean_phases = [PY_mean_phases circ_mean(PY.ripple_phase{un}')];
            
            % significant phase-preference
%             p = circ_otest(PY.ripple_phase{un}); % Hodges-Ajne
            s = sum(PY.ripple_phase{un}<pi/2 & PY.ripple_phase{un}>-pi/2);
            n = s + sum(PY.ripple_phase{un}>pi/2 | PY.ripple_phase{un}<-pi/2);
            p = myBinomTest(s,n,0.5);
            p_PY = [p_PY p];
        end
    end
    for un = 1:numel(IN.ripple_phase)
        if ~isempty(IN.ripple_phase{un})
            IN_mean_phases = [IN_mean_phases circ_mean(IN.ripple_phase{un}')];
            
            % significant phase-preference?
%             p = circ_otest(PY.ripple_phase{un}); % Hodges-Ajne
            s = sum(IN.ripple_phase{un}<pi/2 & IN.ripple_phase{un}>-pi/2);
            n = s + sum(IN.ripple_phase{un}>pi/2 | IN.ripple_phase{un}<-pi/2);
            p = myBinomTest(s,n,0.5);
            p_IN = [p_IN p];
        end
    end
    
    for py = 1:numel(PY.ripple_phase)
        for in = 1:numel(IN.ripple_phase)
            p = circ_wwtest(PY.ripple_phase{py}, IN.ripple_phase{in});
            stats.sig{subj} = [stats.sig{subj}; [p numel(PY.ripple_phase{py}) numel(IN.ripple_phase{in}) wrapTo2Pi(circ_mean(PY.ripple_phase{py}')) wrapTo2Pi(circ_mean(IN.ripple_phase{in}'))] wrapTo2Pi(circ_mean(PY.ripple_phase{py}'))-wrapTo2Pi(circ_mean(IN.ripple_phase{in}'))];
        end
    end
    stats.allSig = [stats.allSig; stats.sig{subj}];
end

[~,~,~,p_FDR] = fdr_bh([p_PY p_IN]);
p_PY_postFDR = p_FDR(1:numel(p_PY));
p_IN_postFDR = p_FDR(numel(p_PY)+1:end);

%%

subject = subjects{1};
for unit_t = 1:numel(unit_types)
    unit_type = unit_types{unit_t};
    load(sprintf('%s/%s/ripple_spike_mat_%s_%s.mat', dataDir, subject, 'DR', unit_type), 'ripple_spike_mat')
    load(sprintf('%s/%s/ripple_spike_mat_LFP_%s_%s.mat',dataDir, subject, 'DR', unit_type), 'ripple_spike_mat_LFP')
        
    raster_img_tmp = zeros(numel(ripple_spike_mat),size(ripple_spike_mat{1},2));
    mean_LFP = NaN(numel(ripple_spike_mat),size(ripple_spike_mat{1},2));
    for un = 1:numel(ripple_spike_mat)
        raster_img_tmp(un,:) = 1000 * mean(ripple_spike_mat{un},1); % instantaneous average spike rate (Hz)
        mean_LFP(un,:) = mean(ripple_spike_mat_LFP{un},1);
    end
    
    % sort raster by largest modulation near t=0
    raster_img = zeros(size(raster_img_tmp));
    rate_vals = mean(raster_img_tmp(:,2001-5:2001+5),2);
    [~,mind] = sort(rate_vals);
    for m = 1:numel(mind)
        raster_img(numel(mind)+1-m,:) = raster_img_tmp(mind(m),:);
    end
    
    figure
    imagesc(raster_img(:,2001-50:2001+50))
    colormap(parula)
    xlabel('time to ripple center (ms)')
    ylabel('unit')
    title([subject ' ' unit_type ' raster'])
    cb = colorbar;
    ylabel(cb, 'mean spike rate (Hz)')
    
    figure
    shadedErrorBar(-50:50,-mean(mean_LFP(:,2001-50:2001+50),1),-nansem(mean_LFP(:,2001-50:2001+50),1));
    xticks(-50:25:50)
    title([subject ' ' unit_type ' raster'])
    ylabel('\muV')
    xlabel('time to ripple center (ms)')
        
end



figure;
hold on
bar(1,100*sum(p_PY_postFDR<0.05)/numel(p_PY_postFDR))
bar(2,100*sum(p_IN_postFDR<0.05)/numel(p_IN_postFDR))
ylabel('percent significant phase-modulation')
ylim([0 100])











