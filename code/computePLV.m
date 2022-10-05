
close all
clear 
clc

addpath(genpath('.')) %entire code directory 

%%

LFP = load(''); %Path to LFP data
rippleStatsPath = ''; %path to rippleStats
load(rippleStatsPath);

subject = 'MG29';
minOverlap = 25; %ms
plvSmoothWin = 10; %ms
win =  500; %ms 
randWin = [-10000 -2000];
nShuff = 200;

RB = zeros(size(LFP));
RBphase = zeros(size(LFP));
[bButter,aButter] = butter(3, [70 100]/(rippleStats.fs/2)); %set up butterworth filter 
for ch = 1:size(LFP,1)
    RB(ch,:) = filtfilt(bButter,aButter,LFP(ch,:));
    RBphase(ch,:) = angle(hilbert(RB(ch,:)));
    fprintf('filterning and calculating phase for channel %i\n', ch)

end

rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
for chRipp = 1:size(rippMask,1) 
    if ~isempty(rippleStats.window{chRipp})
        iS = rippleStats.window{chRipp}(:,1);
        iE = rippleStats.window{chRipp}(:,2);
    
        for ii = 1:length(iE)
            rippMask(chRipp,iS(ii):iE(ii)) = 1;
        end
    end
    
end
%%
PLVstruct = [];
nCh = size(LFP,1);
maxPLV = nan(nCh,nCh);
for chA = 1:nCh
    for chB = 1:nCh
        clc
        fprintf('calculating plv for %i , %i \n', chA, chB)
        
        rippOverlap = rippMask(chA,:) & rippMask(chB,:);
        b = mask2bounds(rippOverlap);
        overlapDuration = b(:,2) - b(:,1);
        b = b(overlapDuration > minOverlap,:);

        coRippleCenters = round(mean(b,2));

        rippA = zeros(length(coRippleCenters), 2*win+1);
        rippB = zeros(length(coRippleCenters), 2*win+1);
        phaseLags = zeros(length(coRippleCenters), 1);

        shufA = zeros(length(coRippleCenters), 2*win+1, nShuff);
        shufB = zeros(length(coRippleCenters), 2*win+1, nShuff);
        for r = 1:length(coRippleCenters)
            cntr = coRippleCenters(r);
            phaseA = RBphase(chA, cntr-win : cntr+win);
            phaseB = RBphase(chB, cntr-win : cntr+win);

            rippA(r,:) = phaseA;
            rippB(r,:) = phaseB;

            phaseAzoom = RBphase(chA, cntr-(win/10) : cntr+(win/10));
            phaseBzoom = RBphase(chB, cntr-(win/10) : cntr+(win/10));
            
            
            phaseLags(r) = circ_mean(angdiff(phaseAzoom, phaseBzoom)');


            for n = 1:nShuff
                randCntr = 0;
                while randCntr - win <= 0 || randCntr + win > rippleStats.recordingLength
                    randCntr = cntr + randi(randWin, 1, 1); % within param.rand_win of each ripple
                end
    
                phaseA = RBphase(chA, randCntr-win : randCntr+win);
                phaseB = RBphase(chB, randCntr-win : randCntr+win);
    
                shufA(r,:,n) = phaseA;
                shufB(r,:,n) = phaseB;
            end
        end

        plv = PLV(rippA', rippB');
        plvSmooth = smoothdata(plv, 'gaussian', plvSmoothWin);
        PLVstruct.plv{chA,chB} = plvSmooth;
        PLVstruct.phaseLags{chA,chB} = plvSmooth;

        for n = 1:nShuff
            plvShuff = PLV(shufA(r,:,n)', shufB(r,:,n)');
            plvSmoothShuff = smoothdata(plvShuff, 'gaussian', plvSmoothWin);

            PLVstruct.plvShuff{chA,chB,n} = plvSmoothShuff;
        end

        maxPLV(chA,chB) = max(plvSmooth);

        
        
    end

           
end

























