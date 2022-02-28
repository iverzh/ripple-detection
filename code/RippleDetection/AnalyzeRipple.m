 function [rippleSleepSet, rejectVec] = AnalyzeRipple(sleepdata, rippleband, win, rippleWindow, IIScheck, spikeMask, nan_edge_mask, ...
                                                      IISflag, ch, goodRipples, RejectParams, fs, chan_labels, shift, rejectFlag)

% notch filter data to remove line noise (60 Hz + harmonics)
% for nf = 60:60:fs/2 %240
%     Wo = nf/(fs/2);
%     BW = Wo/35;
%     [b,a] = iirnotch(Wo, BW);
%     sleepdata = filtfilt(b,a,sleepdataOrig);
% end

sleepdata = sleepdata(ch,:);

% bandpass at ripple band (70-100 Hz)
[b,a] = butter(3,rippleband/(fs/2));
ripple.rippleband = filtfilt(b,a,sleepdata); % zero phase

ind = find(ripple.rippleband>500);
for ii = ind
    if ii < fs  
       ripple.rippleband(1:ii+fs) = 0;
       sleepdata(1:ii+fs)         = 0;
    elseif ii+fs > length(ripple.rippleband)
       ripple.rippleband(ii-fs:end) = 0;
       sleepdata(ii-fs:end)         = 0;
    else
       ripple.rippleband(ii-fs:ii+fs) = 0;
       sleepdata(ii-fs:ii+fs)         = 0;
    end
end

% bandpass at low frequency band (25-50 Hz) 
[b,a] = butter(3,100/(fs/2),'high');
ripple.sharpspike = filtfilt(b,a,sleepdata); % zero phase


% bandpass at high frequency band (200- Hz) 
% [b,a] = butter(3,200/(fs/2), 'high');
[b,a] = butter(3,200/(fs/2), 'high');
ripple.highfreqband = filtfilt(b,a,sleepdata); % zero phase

% bandpass at middle frequency band (20- Hz) (too detect large epileptic events in LFP) 
% [b,a] = butter(3,150/(fs/2), 'high');
[b,a] = butter(3,150/(fs/2), 'high');
ripple.HGband = filtfilt(b,a,sleepdata); % zero phase

% bandpass at low frequency band (25-50 Hz) 
[b,a] = butter(3,[100 200]/(fs/2));
% [b,a] = butter(3,100/(fs/2),'high');
ripple.band100200 = filtfilt(b,a,sleepdata); % zero phase


a = sprintf('notch and bandpass filter channel %i\n', ch);
fprintf(a)

zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

nanData = abs(hilbert(ripple.rippleband));
nanData(nan_edge_mask) = NaN;
ripple.RBzscore = zscore_xnan(nanData); %zscore of env of hilbert amplitudes

nanData = abs(hilbert(ripple.highfreqband));
nanData(nan_edge_mask) = NaN;
ripple.HFzscore = zscore_xnan(nanData);

nanData = abs(hilbert(ripple.sharpspike));
nanData(nan_edge_mask) = NaN;
ripple.sharpzscore = zscore_xnan(nanData);

nanData = abs(hilbert(ripple.HGband));
nanData(nan_edge_mask) = NaN;
ripple.HGzscore = zscore_xnan(nanData);

nanData = abs(hilbert(ripple.band100200));
nanData(nan_edge_mask) = NaN;
ripple.zscore100200 = zscore_xnan(nanData);

nanData = sleepdata;
nanData(nan_edge_mask) = NaN;
ripple.diffLFPzscore  = abs(zscore_xnan(diff(nanData))); %zscore of LFP numerical derivative 

ripple.diffLFP  = abs(diff(sleepdata)); %zscore of LFP numerical derivative 

ripple.LFP = sleepdata; 

ripple.ind = goodRipples;
           
ripple.ind = ripple.ind(ripple.ind- (win+rippleWindow) > 0 & ripple.ind+ (win+rippleWindow) < size(ripple.LFP,2));
ripple = recenterRipples(ripple,rippleWindow, win);
ripple.ind = ripple.ind(ripple.ind-win > 0 & ripple.ind+win < size(ripple.LFP,2));

mergedLocs = mergeRipples(ripple,0.025,fs);
ripple.ind = mergedLocs;



ripple.events.raw               = NaN(numel(ripple.ind), win*2 + 1);
ripple.events.rippleband        = NaN(numel(ripple.ind), win*2 + 1);
ripple.events.highFreqBand      = NaN(numel(ripple.ind), win*2 + 1);
ripple.events.band100200        = NaN(numel(ripple.ind), win*2 + 1);
% ripple.events.midFreqBand       = NaN(numel(ripple.ind), win*2 + 1);

ripple.events.RBzscore          = NaN(numel(ripple.ind), win*2 + 1);
ripple.events.HFzscore          = NaN(numel(ripple.ind), win*2 + 1);
ripple.events.sharpzscore       = NaN(numel(ripple.ind), win*2 + 1);
ripple.events.HGzscore          = NaN(numel(ripple.ind), win*2 + 1);
ripple.events.zscore100200      = NaN(numel(ripple.ind), win*2 + 1);
ripple.events.diffLFP           = NaN(numel(ripple.ind), win*2 + 1);  
ripple.events.diffLFPzscore     = NaN(numel(ripple.ind), win*2 + 1);  

ripple.events.rippleAmp         = NaN(numel(ripple.ind),1);

fprintf(['extracting ripples and ripple stats from channel : ',num2str(chan_labels{ch}), '\n'])      


for rip_num = 1:numel(ripple.ind) % loop through ripple events

    ripple.events.raw(rip_num,:)          = ripple.LFP(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);
    ripple.events.rippleband(rip_num,:)   = ripple.rippleband(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);
    ripple.events.highFreqBand(rip_num,:) = ripple.highfreqband(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);
    ripple.events.band100200(rip_num,:)   = ripple.band100200(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);
%                 ripple.events.midFreqBand(rip_num,:)  = ripple.band100200(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);

    ripple.events.HFzscore(rip_num,:)     = ripple.HFzscore(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);   
    ripple.events.sharpzscore(rip_num,:)  = ripple.sharpzscore(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);   
    ripple.events.RBzscore(rip_num,:)     = ripple.RBzscore(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);   
    ripple.events.HGzscore(rip_num,:)     = ripple.HGzscore(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);   
    ripple.events.zscore100200(rip_num,:) = ripple.zscore100200(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);   
    ripple.events.diffLFP(rip_num,:)      = ripple.diffLFP(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);   
    ripple.events.diffLFPzscore(rip_num,:)= ripple.diffLFPzscore(ripple.ind(rip_num)-win:ripple.ind(rip_num)+win);   

    % compute ripple stats
    % 20Hz Highpass range
%     midfreqRange = range(ripple.events.midFreqBand(rip_num,:));
%     ripple.events.midFreqRange(rip_num)        = midfreqRange;

    center = size(ripple.events.raw,2)/2;
    center = round(center);
    rippleRange = center-rippleWindow:center+rippleWindow;

    %High Frequency Hilbert Amp
%     signalhighF =ripple.events.highFreqBand(rip_num,rippleRange);
%     ripple.events.highFreqAmp(rip_num)       = max(abs(hilbert(signalhighF)));

    % Low Frequency Ripple Amp
    signallowF =ripple.events.band100200(rip_num,rippleRange);
    ripple.events.lowFreqAmp(rip_num)        = max(abs(hilbert(signallowF)));

    % Ripple Band Amp
    signalRipple =ripple.events.rippleband(rip_num,center);
    ripple.events.rippleAmp(rip_num)         = max(abs(hilbert(signalRipple)));


end

if rejectFlag
    [tempStruct, rejectVec] = rejectRipples(ripple, RejectParams, rippleWindow, ch, fs, IISflag, spikeMask, IIScheck);
else
    tempStruct = ripple.events;
    tempStruct.goodRipples = ripple.ind;
    rejectVec = [];
end

fn = fieldnames(tempStruct);
for f = 1:numel(fn)
    if isfield(ripple.events, fn{f})
        ripple.events.(fn{f}) = []; %delete for memory
    end
    rippleSleepSet.(fn{f}) = tempStruct.(fn{f});
end


rippleSleepSet.goodRipplesConcat = rippleSleepSet.goodRipples + shift; %shift indicies to account for previous sleep sets


return




