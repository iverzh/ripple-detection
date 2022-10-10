function rippleStatsFinal =  computeRippleStats(events, RejectParams, fs)


smthKrnl    = RejectParams.bndParams.smthKrnl;
srchWinSz   = RejectParams.bndParams.srchWinSz;
scoreThresh = RejectParams.bndParams.scoreThresh;


if iscell(events.raw)
    nChan = size(events.raw,2); 
else
    nChan = 1;
end

% pre-allocate
plotStats.centeredInd    = cell(1,nChan);
plotStats.duration       = cell(1,nChan);
plotStats.HGz            = cell(1,nChan);
plotStats.InterRipPeriod = cell(1,nChan);
plotStats.oscFreq        = cell(1,nChan);
plotStats.rippleAmp      = cell(1,nChan);
plotStats.window         = cell(1,nChan);

for ch = 1:nChan
    nRip = length(events.goodRipples);
    centered = ones(nRip,1);
    if nRip > 0
        for rip = 1:nRip

            %Ripple Amplitude
            plotStats.rippleAmp{ch}(rip) = events.rippleAmp(rip);

            %Ripple zero cross
            center = round(size(events.raw(rip, :),2)/2);
            rippleBand = events.rippleband(rip, :);
            RBzscore = events.RBzscore(rip,:);
            RBzscoreSmooth = smoothdata(RBzscore,2,'gaussian',smthKrnl);

            HGzscore = events.HGzscore(rip,:);

%             score = Inf;
%             winStart = 0;
%             while score > 0.75 && winStart < 200
%                 score = abs(RBzscore(center-winStart));
%                 winStart = winStart + 1;
%             end
% 
%             score = Inf;
%             winEnd = 0;
%             while score > 0.75 && winEnd < 200 
%                 score = abs(RBzscore(center+winEnd));
%                 winEnd = winEnd + 1;
%             end

            score = Inf;
            winStart = 0;
            while score > scoreThresh && winStart < fs*2%200
                score = RBzscoreSmooth(center-winStart);
                winStart = winStart + 1;
            end
            range=(-round(smthKrnl*(srchWinSz/2)):round(smthKrnl*(srchWinSz/2)))+center-winStart;
            range(range>center)=[];
            range(range<1)=[];
            %[~,idx] = min(RBzscore(range));
            if all(isnan(RBzscore(range)))
                winStart= (center-range(end)) - find(~isnan(RBzscore(range(end):center)),1,'first');
            elseif winStart~=fs*2
                idx=[];
                bump=0;
                while isempty(idx)
                    idx=find(RBzscore(range)<scoreThresh+bump,1,'last');
                    bump=bump+0.1;
                end
                winStart = winStart - (idx - round(smthKrnl*(srchWinSz/2)));
            end
            
            score = Inf;
            winEnd = 0;
            while score > scoreThresh && winEnd < fs*2 %200
                score = RBzscoreSmooth(center+winEnd);
                winEnd = winEnd + 1;
            end
            range=(-round(smthKrnl*(srchWinSz/2)):round(smthKrnl*(srchWinSz/2)))+center+winEnd;
            nTrim=sum(range<center);
            range(range<center)=[];
            range(range>length(RBzscore))=[];
            %[~,idx]=min(RBzscore(range));
            if all(isnan(RBzscore(range)))
                winEnd=find(~isnan(RBzscore(center:range(1))),1,'last');
            elseif winEnd~=fs*2
                idx=[];
                bump=0;
                while isempty(idx)
                    idx=find(RBzscore(range)<scoreThresh+bump,1,'first');
                    bump=bump+0.1;
                end
                winEnd = winEnd + (idx - (round(smthKrnl*(srchWinSz/2)) - nTrim));
            end

            plotStats.duration{ch}(rip) = (winStart + winEnd) * (1/fs) * 1000; %ms

            loc = events.goodRipplesConcat(rip);
            plotStats.window{ch}(rip,:) = [loc-winStart, loc+winEnd];
            plotStats.HGz{ch}(rip) = mean(HGzscore(center-winStart:center+winEnd));

            window = center-winStart:center+winEnd;
            winData = rippleBand(window);
            phases = angle(hilbert((winData))); % / pi;
            
            zcs = find(winData(1:end-1).*winData(2:end) < 0); %zero crossings
            
            if ~isempty(zcs)
                % count additional oscillation that extends beyond zero
                % crossings
                unWrapStart = unwrap(phases(1:zcs(1)));
                add1 =  ( unWrapStart(end) - unWrapStart(1) ) / pi;

                unWrapEnd = unwrap(phases(zcs(end):end));
                add2 =  ( unWrapEnd(end) - unWrapEnd(1) ) / pi;


                nCycles  = (length(zcs) + add1 + add2) /2;



                freq = round( nCycles / (length(window)/fs)*10)/10;
            else
                freq = 0;
            end
        
        
     

            plotStats.oscFreq{ch}(rip) = freq;

            if plotStats.duration{ch}(rip) < 3*(1/100 * fs) %80 Hz is putative ripple frequency 
                keep(rip) = 0;
            end

            
        end    
    end
    
    plotStats.InterRipPeriod{ch} = diff(events.goodRipplesConcat) / fs ;

% %     centered(centered==0)=1;
%     centered = logical(centered);
%     if length(centered) > 1
%         rippleStatsFinal.InterRipPeriod{ch} = diff(events.goodRipplesConcat{ch}(centered)) / fs ;
%         rippleStatsFinal.duration{ch} = plotStats.duration{ch}(centered);
%         rippleStatsFinal.oscFreq{ch} = plotStats.oscFreq{ch}(centered);
%         rippleStatsFinal.rippleAmp{ch} = plotStats.rippleAmp{ch}(centered);
%         rippleStatsFinal.centeredInd{ch} = centered;
%         rippleStatsFinal.window{ch} = plotStats.window{ch};
%         rippleStatsFinal.HGz{ch} = plotStats.HGz{ch};
%     end    


end

rippleStatsFinal = plotStats;



return