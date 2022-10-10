
function [ii, data, fs] = DetectHighFreqEvents(BroadbandData,bandpassRange,RejectParams,fs,nan_edge_mask,min_cycles,zs_thresh)

zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

if diff(bandpassRange) < 150
    data = BroadbandData;
    fsOrig = fs;
else
    data = resample(BroadbandData, 1e4, fs);
    fsOrig = fs;
    fs = 1e4;
end

if any(isnan(data))
    nanMask = isnan(data); 
    data(nanMask) = 0; 
else
    nanMask = false(1,length(data));
end

[b,a] = butter(3,bandpassRange/(fs/2));
bandpass = filtfilt(b,a,data); % zero phase

nanData = abs(hilbert(bandpass));
nanData(nan_edge_mask) = NaN;
nanData(nanMask) = NaN;
bandpass_zscore = zscore_xnan(nanData); % zscore of env of hilbert amplitudes


HFregions = mask2bounds(bandpass_zscore > zs_thresh); %get bounds of all high frequency events

del = zeros(size(HFregions,1),1);
% remove events that are less than 3 cycles
for iR = 1:size(HFregions,1)
    winData = bandpass(HFregions(iR,1):HFregions(iR,2)); 
    zcs = find(winData(1:end-1).*winData(2:end) < 0); %zero crossings
    
    if length(zcs)/2 < min_cycles; del(iR) = 1; end

end

HFregions(logical(del), :) = [];

HFdurations = HFregions(:,2) - HFregions(:,1);
HFregions(HFdurations < RejectParams.minDuration(1)*fs, :) = []; %delete events that are too short 

% ii = round(mean(HFregions,2) / (fs/fsOrig))';
ii = round(mean(HFregions,2))';

return 


















