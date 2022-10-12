function [spike_mask, spike_regs] = detectSpikes_IAV(dat, fs, hge, hf_factor, spike_factor, spike_score_cutoff, spike_pad, hf_shift)
% v2 edit to adjust for sampling rate
nanMask = isnan(dat(:,1));
dat(nanMask,:) = [];
% CHARLIE: I used the following parameters to good effect when I tested
% this on your cortical high gamma events:
if nargin < 8
    hf_shift = round(0.04*fs); %samples
end
if nargin < 7
    spike_pad = 0*fs; %pad the spike on each side in the returned mask
end
if nargin < 6
    spike_score_cutoff = 130; %cutoff for the weighted sum of hf and spike scores for detecting a spike
end
if nargin < 5
    spike_factor = 5/(fs/500); %weight for the spike score (corrected for convolution scaling with num samples)
end
if nargin < 4
    hf_factor = 13; %weight for the high frequency score
end
if nargin < 3 || isempty(hge)
    disp('high gamma envelope...')
    if fs/2 > 190
        [b,a] = butter(4, [70 190]/(fs/2), 'bandpass');
    else
        [b,a] = butter(4, 70/(fs/2), 'high');
    end
    
        hge = abs(hilbert(filtfilt(b,a,dat)));
    hge = uint16(hge*100);
    %multiplied by 100 to preserve acceptable resolution once converted to uint16;
end %change last argument to zero to reduce memory overhead

num_channels = size(dat,2);

%% high frequency score

disp('smoothing high gamma envelope...')
tic
hf_score = uint16(zeros(size(dat)));
for n = 1:num_channels
    if num_channels > 1
        disp(n/num_channels)
    end
    hf_score(:,n) = uint16(smooth(single(hge(:,n)),20));
end
if hf_shift > 0
    hf_score = [zeros(hf_shift,num_channels); hf_score(1:end-hf_shift,:)];
elseif hf_shift < 0
    hf_score = [hf_score(-hf_shift:end,:); zeros(-hf_shift,num_channels)];
end
hf_med = double(median(hf_score));
toc

% Spike convolution

load('/home/jgarret/ArtifactDetection/wavelets_detrended.mat', 'spike_wavelet2')
if fs ~= 500
    spike_wavelet2 = resample(spike_wavelet2,fs,500);
end

disp('computing spike template convolution...')
tic
spike_score = uint16(zeros(size(dat)));
for n = 1:num_channels
    spike_score(:,n) = uint16(abs(conv(dat(:,n), flip(spike_wavelet2), 'same'))*100);
end
spike_score = padToLength(spike_score, length(dat));
spike_score_med = double(median(spike_score));
toc

% parameters

hf_weight = hf_factor./double(hf_med);
spike_weight = spike_factor./double(spike_score_med);

% find spikes

spike_regs = cell(num_channels,1);
for n = 1:num_channels
    disp(n/num_channels)
    spike_regs{n} = mask2bounds(spike_score(:,n)*spike_weight(n) + hf_score(:,n)*hf_weight(n) > spike_score_cutoff);
end

% spike mask

spike_mask_clipped = false(size(dat));
for n = 1:num_channels
    spike_mask_clipped(:,n) = bounds2mask(spike_regs{n}, length(dat), spike_pad);
end

ii = find(~nanMask);
spike_mask = false(size(dat,2),length(nanMask));
 for n = 1:num_channels
    ii_n = ii(spike_mask_clipped(:,n));
    spike_mask(:,ii_n) = true;

 end



disp('done')

end


%% SUBFUNCTION

function mask = bounds2mask(bounds, mask_len, num_pad)
%BOUNDS2MASK takes an n x 2 matrix representing the bounds of n regions in
%a 1d array and returns a logical mask of length mask_len representing
%those regions, with each region optionally padded on both sides by num_pad
%elements
if nargin < 3
    num_pad = 0;
end

mask = false(mask_len,1);
if ~isempty(bounds)
    for n = 1:size(bounds,1)
        curr_bounds = bounds(n,:);
        curr_bounds(1) = curr_bounds(1) - num_pad;
        curr_bounds(2) = curr_bounds(2) + num_pad;
        curr_bounds(curr_bounds<1) = 1;
        curr_bounds(curr_bounds>mask_len) = mask_len;
        if(curr_bounds(1)<=curr_bounds(2))
            mask(curr_bounds(1):curr_bounds(2)) = true;
        end
    end
end
end

