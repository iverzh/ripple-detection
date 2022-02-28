function [ SOmat, summary ] = detectSO( data, sfreq, clust_mask, per)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect Slow Oscillations
% Taken from Jason Trees detect_slow_waves.m 8.5.15
% Func implementation by CEG 6.16.16

% INPUT_ARGS:
% data= (# channel x time) double
% sfreq = sampling frequency;
% clust_mask = binary mask length of data (enter '1' for samples of interest). 
% per = Top and bottom percent of peaks to include (e.g. per = 40);
% NOTE: This script assumes down states are depth-negative. 


% OUTPUT_ARGS:
% SOmat = (# channel* # SO x 3 col) where
%         col1-3 = sample of peak/trough, amplitude, channel #
% summary = (# channel x 3 col) where
%         col1-3 = # Up and Down states, # bad indices (e.g. gaps or
%         sharp), # events excluded because close to bad indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if nargin < 3
    error('Too few inputs');
elseif nargin < 4
    out_path = [];
end

rej_thresh = 15; %number of std above mean of abs(data) used for artifact threshhold
%per = 40; %outer bound percentage of peaks to keep
filter_range = [0.1 4];
lat_range = [0.25 3]; %latency of zero consecutive zero crossings (in seconds)

all_slow_peaks = [];
recap_mat = [];
for ch = 1:size(data,1)
    
    ch_data = data(ch,:);
    
    % Find diff artifact indices (due to cleveland gaps) 
    art_ind=[];
    absdat = abs(diff(ch_data));
    mdat = mean(absdat);
    sdat = std(absdat);
    art_ind = find(absdat>mdat+(rej_thresh*sdat));
    
%     %Sanity plot
%     figure('Visible','on');
%     plot(absdat)
%     hold on
%     plot([0,length(absdat)],[mdat+(rej_thresh*sdat),mdat+(rej_thresh*sdat)],'m:')
%     plot([0,length(absdat)],[mdat+(rej_thresh*sdat),mdat+((rej_thresh+1)*sdat)],'m:')
%     
%     figure('Visible','on');
%     test_out = find(absdat>9 & absdat<20);
%     s = art_ind(3);
%     plot(ch_data(1,s-(10*sfreq):s+(10*sfreq)));

    % High pass to remove low freq component
    [b,a] = butter(4, 20/(sfreq/2),'high');
    art_filt = filtfilt(b,a,double(ch_data));
    absfilt = abs(art_filt);
    
    % Calculate artifactual deviations of the amplitude
    % NOTE: This is kind of redundant on top of the first pass, but okay...
    mfdat = mean(absfilt);
    sfdat = std(absfilt);
    art_filt_ind = find(absfilt > mfdat+((rej_thresh)*sfdat));
%     figure('Visible','on');
%     s = art_filt_ind(5);
%     plot(ch_data(1,s-(10*sfreq):s+(10*sfreq)));
    
    art_ind = [art_ind art_filt_ind];
    
    % Filter
    if length(filter_range) == 2
        fprintf('Filtering from %.1d to %d Hz...\n',filter_range);
        %lowpass
        [b,a] = butter(4,filter_range(2)/(sfreq/2),'low');
        tmp_filt1 = filtfilt(b,a,double(ch_data));
        %highpass
        [b,a] = butter(4,filter_range(1)/(sfreq/2),'high');
        tmp_filt = filtfilt(b,a,double(tmp_filt1));
    else
        fprintf('Lowpass at %d Hz...\n',filter_range);
        %lowpass
        [b,a] = butter(4,filter_range/(sfreq/2),'low');
        tmp_filt = filtfilt(b,a,double(ch_data));
    end;
    clear tmp_filt1 
%     figure('Visible','on');
%     plot(ch_data(1,1*sfreq:100*sfreq));
%     hold on
%     plot(tmp_filt(1,1*sfreq:100*sfreq),'r')
    
    % Find zero crossings ind
    indx = find(tmp_filt(1:end-1).*tmp_filt(2:end)<0);
    
    % Find indices within specified range
    tmp_slow_ind = find(diff(indx)>=lat_range(1)*sfreq & diff(indx)<=lat_range(2)*sfreq);
    slow_ind_bounds = [indx(tmp_slow_ind)',indx(tmp_slow_ind+1)'];
    
    % Find max peak in range
    slow_ind_peaks = [];
    for sp = 1:length(slow_ind_bounds);
        rel_ind = find(abs(tmp_filt(slow_ind_bounds(sp,1):slow_ind_bounds(sp,2))) == max(abs(tmp_filt(slow_ind_bounds(sp,1):slow_ind_bounds(sp,2)))));
        slow_ind_peaks = [slow_ind_peaks, slow_ind_bounds(sp,1)+rel_ind-1];
    end
    
    % Eliminate points not in mask
    if length(clust_mask)>0
        elim_peaks = find(clust_mask(slow_ind_peaks)==0);
        slow_ind_peaks(elim_peaks) = [];
    end
    
    % Eliminate points too close to artifact indices (2s)
    art_peaks = [];
    for jj = 1:length(slow_ind_peaks)
        tmp_art_ind = abs(art_ind - slow_ind_peaks(jj));
        tmp_art_val = tmp_art_ind(tmp_art_ind == min(tmp_art_ind));
        if tmp_art_val <= 2*sfreq
            art_peaks = [art_peaks, jj];
        end;
    end;
    slow_ind_peaks(art_peaks) = [];
    
    % Select out only outer limits of peaks
    peak_val = tmp_filt(slow_ind_peaks);
    pos_ind = find(peak_val>0);
    neg_ind = find(peak_val<0);
    keep_pos = find(peak_val(pos_ind)> prctile(peak_val(pos_ind),100-per));
    keep_neg = find(peak_val(neg_ind)< prctile(peak_val(neg_ind),per));
    slow_ind_peaks = slow_ind_peaks(sort([pos_ind(keep_pos),neg_ind(keep_neg)]));
    
    % Concatenate the peak indices with channels
    all_slow_peaks = [all_slow_peaks; [slow_ind_peaks',tmp_filt(slow_ind_peaks)', ch*ones(size(slow_ind_peaks))']];
    recap_mat = [recap_mat; [length(slow_ind_peaks) length(art_ind) length(art_peaks)]];
    close all
    
%     % Check some
%     time_to_down = -4:1/sfreq:4;
%     ncheck = 5;
%     DS = all_slow_peaks(all_slow_peaks(:,2)<0,1);
%     figure('Visible','on');
%     for pp = 1:ncheck
%         down = DS(randperm(length(DS)));
%         jitter = 900*pp;
%         
%         plot(time_to_down,ch_data(1,down-4*sfreq:down+4*sfreq)*8+jitter);
%         hold on
%     end
%     hold off
%     
end;

SOmat = all_slow_peaks; summary = recap_mat;
summary

et = toc;
display(['Total elapsed time is ' num2str(et/60) ' minutes.']);
end

