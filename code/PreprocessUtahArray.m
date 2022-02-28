% PreprocessUtahArray.m 
% Script used to downsample and remove unit action potentials from broadband
% Utah Array data in the following manuscript:
% 
% Widespread ripples synchronize human cortical activity during sleep, waking, 
% and memory recall. CW Dickey, IA Verzhbinsky, X Jiang, BQ Rosen, S Kajfez, 
% B Stedelin, JJ Shih, S Ben-Haim, AM Raslan, EN Eskander, J Gonzalez-Martinez,
% SS Cash, E Halgren
% 
% Ilya A. Verzhbinsky, Halgren Lab, 02.24.2022

clc
clear 
close all

addpath(genpath('./'))
%% Inputs

subject = 'U1'; %name of subject

rawDataPath = './'; %path to the directory of the raw .ns5 utah array files. 
spikePath = './../data/units.mat'; %path to spike time structure. See ./data/units.mat for example format. 
outPath = './../data/UA/'; %path to output folder for downsampled Utah Array data. 
if ~isfolder(outPath); mkdir(outPath); end

flist = dir(fullfile(rawDataPath,'*.ns5'));
flist = {flist.name};

fsDownsample = 1000; %target downsample rate in Hz;
fs = NS5.MetaTags.SamplingFreq;
startTime = 15; % time in samples before unit peak to begin computing the unit template
endTime = 50; % time in samples after unit peak to end computing the unit template
recenterWin = 60;

nSessions = length(flist);
%% Load Raw Data
rawData = [];
for ns = 1:nSessions
    clc
    if ismember(sessionID(ns),sessions)
        fprintf('loading file %s (%i out of %i)\n',flist{ns},ns,nSessions)

        NSX_open(fullfile(rawDataPath,flist{ns}),'read','report')

        rawData = [rawData, NS5.Data];
    end
end

clear NS5.Data

AR = false; %Average reference utah data.
%%
if strcmp(NS5.MetaTags.FileTypeID,'NEURALSG')
    
elseif strcmp(NS5.MetaTags.FileTypeID,'NEURALCD')
    dig_min = double(NS5.ElectrodesInfo(1).MinDigiValue);
    dig_max = double(NS5.ElectrodesInfo(1).MaxDigiValue);
    phy_min = double(NS5.ElectrodesInfo(1).MinAnalogValue);
    phy_max = double(NS5.ElectrodesInfo(1).MaxAnalogValue);

    % Convert to phyical units
    rawData = double(rawData);
    rawData = (rawData-dig_min)./(dig_max-dig_min);
    rawData = rawData.*(phy_max-phy_min)+phy_min;  
%     
%         tmp = (tmp-dig_min)./(dig_max-dig_min);
% %     tmp = tmp.*(phy_max-phy_min)+phy_min;  
end

dataLength = length(rawData) / NS5.MetaTags.SamplingFreq / 60/ 60; 
minutes = round(dataLength * 60);
fprintf('%f hours (%i minutes) of raw data loaded\n',dataLength, minutes)

%% Load spikes
close all

load(spikePath)

dataUnitsRemoved = rawData;
clear rawData

%% Compute spike template and subtract our unit waveforms
spikeTemplates = zeros(size(units,1), startTime+endTime+1);
for ch = 1:size(units,1)
    
    switch subject 
        case 'MG67_'
            fileName = units{ch,1};
            chanLabel = str2double(fileName(18:19));
            if isnan(chanLabel); chanLabel = str2double(fileName(18)); end
            
        case {'MG29','MG49','MG67'}
            chanLabel = units{ch,1};
    end
% 
%     % Convert to phyical units
    dataChanRemUnits = double(dataUnitsRemoved(chanLabel,:));
    spikeTimes = units{ch,2} * fs;
    spikeTimes(spikeTimes > (length(dataChanRemUnits) - recenterWin)) = [];
    spikeTimes(spikeTimes < recenterWin) = [];
    nSpk  = length(spikeTimes);
    
    if nSpk > 0
        
        spikeMat = zeros(nSpk,startTime+endTime+1);
        for spk = 1:nSpk
            times = round(spikeTimes(spk)-recenterWin:spikeTimes(spk)+recenterWin); %times for recentering 

            spikeOrig = dataChanRemUnits(times);
            shift = find(spikeOrig == min(spikeOrig), 1, 'first'); %
            shift = shift - round((2*recenterWin+1)/2);
        %     
            spikeTimes(spk) = spikeTimes(spk) + shift;

            times = round(spikeTimes(spk)-startTime:spikeTimes(spk)+endTime);
            spikeRecenter = dataChanRemUnits(times);
            spike = spikeRecenter - spikeRecenter(1);

            spikeMat(spk,:) = spike;
        end

        remSpikes = mean(spikeMat);


        for spk = 1:nSpk
            times = round(spikeTimes(spk)-startTime:spikeTimes(spk)+endTime);

            spikeOrig = dataChanRemUnits(times);
            spike = spikeOrig;
            spike = spike  -  remSpikes; %subtract template
            spikeMat(spk,:) = spike;
            dataChanRemUnits(times) = spike;
        end
            
        

        dataUnitsRemoved(chanLabel,:) = dataChanRemUnits;
        spikeTemplates(ch,:) = remSpikes;
        fprintf('%i units cut out\n', chanLabel)
    end
end



%% Down sample raw data


gMsk = ones(1,size(dataUnitsRemoved,1));

LFP = [];
count = 1;
for ch = 1:size(dataUnitsRemoved,1)
    tmp = double(dataUnitsRemoved(ch,:));
    if ismember(ch,gMsk)
        LFP(count,:) = resample(tmp,fsDownsample,fs);
        
        % notch filter data to remove line noise (60 Hz + harmonics)
        for nf = 60:60:fsDownsample/2 %240
            Wo = nf/(fsDownsample/2);
            BW = Wo/35;
            [b,a] = iirnotch(Wo, BW);
            LFP(count,:) = filtfilt(b,a,LFP(count,:));
        end

        fprintf([num2str(ch),' done \n'])
        count = count + 1;
    end
end


if AR %compute average reference 
    ref = mean(LFP);
    LFP = LFP - ref;
end

save(fullfile(outPath,[subject, '_downsample.mat']))







