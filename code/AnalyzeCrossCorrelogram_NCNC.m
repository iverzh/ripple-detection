% AnalyzeCrossCorrelogram.m 
% Script used to produce all neocortical-neocrotical ripple cross correlogram
% statistics and plots in the following manuscript:
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



addpath(genpath('.')) %entire code directory 


%% Inputs

subj_list_full = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
                  'S11','S12','S13','S14','S15','S16','S17'};


runList        = 1:17; %list of subjects you want to run
state          = 'wake'; %wake or sleep

PRTHfolderNC   = [ './../data/CrossCorrelogram/',state,'/NCNC'];
badChanDir     = './../data/bad_chan';
exportDir      =  './../Figures';
if ~isfolder(exportDir); mkdir(exportDir); end


FDRwin        = 500; % +- window (in ms) to compute FDR significance
statsBinWidth = 25;  %histogram bin width (in ms) for computing all stats (significance, sidedness)
plotBinWidth  = 25; %histogram bin width (in ms) for plotting
fs            = 1000; %sample rate

%% load peri-ripple time histogram data

NCp = [];
NCpFDR = [];
NCa = [];
NCc = [];
NCn = [] ; 

NCa_plot = [];
NCn_plot = [];

subjNames = {};
subjID = [];

c = 0;
nRip = 0;
for subj = runList %length(subj_list_full)
    subject = subj_list_full{subj};

    NC = load(fullfile(PRTHfolderNC,[subject,'_NCNC_PRTH.mat']));
    
    % compute histogram parameters
    nIter = size(NC.subjPRTH.nullPRTH{1,2},1);
    win = ceil(size(NC.subjPRTH.nullPRTH{1,2},2)/2);
    edges = -win+(statsBinWidth/2):statsBinWidth:win;
    edgesPlot = -win+(plotBinWidth/2):plotBinWidth:win;
    times = -win:win;
    nBins = 2*win/statsBinWidth - 1;
    nBinsPlot = 2*win/plotBinWidth - 1;
    center = ceil(nBins/2);


    load(fullfile(badChanDir,[subject,'_bad_chan.mat']))
    
    
    for ch1 = 1:length(NC.subjPRTH.chanLabels)
        chanLabel1 = NC.subjPRTH.chanLabels{ch1};
        
        for ch2 = 1:length(NC.subjPRTH.chanLabels)
            chanLabel2 = NC.subjPRTH.chanLabels{ch2};

            if ~contains(chanLabel1,bad_chan) && ~contains(chanLabel2,bad_chan) && ch1 ~= ch2

                countsTemp = smoothdata(NC.subjPRTH.eventPRTH{ch1,ch2}, 'gaussian',250);
                nullTemp   = smoothdata(NC.subjPRTH.nullPRTH{ch1,ch2},2, 'gaussian',250);
                counts = zeros(1,nBins);
                null = zeros(nIter,nBins);

                for e = 1:length(edges)-1
                   ii = times > edges(e) & times < edges(e+1);
                   counts(e) = sum(countsTemp(ii));
                   null(:,e) = sum(nullTemp(:,ii),2);
                end
                
                countsTemp = smoothdata(NC.subjPRTH.eventPRTH{ch1,ch2}, 'gaussian',50);
                nullTemp   = smoothdata(NC.subjPRTH.nullPRTH{ch1,ch2}, 2, 'gaussian',50);
                countsOrig = zeros(1,nBinsPlot);
                nullOrig = zeros(nIter,nBinsPlot);

                for e = 1:length(edgesPlot)-1
                   ii = times > edgesPlot(e) & times < edgesPlot(e+1);
                   countsOrig(e) = sum(countsTemp(ii));
                   nullOrig(:,e) = sum(nullTemp(1:nIter,ii),2);
                end

                countsPlot = countsOrig; %smoothdata(countsOrig, 'gaussian',5);
                nullPlot   = nullOrig; %smoothdata(nullOrig, 2, 'gaussian',5);



                n =  NC.subjPRTH.nRipples{ch1,ch2};
                NCp = [NCp, findPercentile(null,counts,'pval')];
                NCpFDR = [NCpFDR, findPercentile(null(:,center-(FDRwin/statsBinWidth):center+(FDRwin/statsBinWidth)),counts(center-(FDRwin/statsBinWidth):center+(FDRwin/statsBinWidth)),'pval')];

                NCa = [NCa, counts/(n*statsBinWidth)*fs]; 
                NCa_plot = [NCa_plot, countsPlot/(n*plotBinWidth)*fs]; 


                NCc = [NCc, NC.subjPRTH.eventPRTH{ch1,ch2}]; 

                NCn = [NCn, null/(n*statsBinWidth)*fs]; 
                NCn_plot = [NCn_plot, nullPlot/(n*plotBinWidth)*fs]; 

                c = c + 1;
                subjNames{c} = subject;
                subjID(c) = subj;

            
            end
        end
        
%         NCp = [NCp, NC.subjPRTH.(mode).nullPRTH{ch}]; 
    end
    
    fprintf([subject,' done. \n'])

end


%% Compute FDR
FDRbins = 2*FDRwin/statsBinWidth+1;

[hNC, ~, ~, adjNCp]  = fdr_bh(NCpFDR, 0.05, 'pdep','yes');


adjNCp = reshape(adjNCp, [length(adjNCp)/FDRbins, FDRbins]);
hNC    = reshape(hNC,    [FDRbins, length(hNC)/FDRbins]);

NCa = reshape(NCa, [nBins, length(NCa)/nBins]);
NCa_plot = reshape(NCa_plot, [nBinsPlot, length(NCa_plot)/nBinsPlot]);

NCc = reshape(NCc, [2*win-1, length(NCc)/(2*win-1)]);

NCn = reshape(NCn, [nIter, nBins, length(NCn)/nBins]);
NCn_plot = reshape(NCn_plot, [nIter, nBinsPlot, length(NCn_plot)/nBinsPlot]);

adjNCp = adjNCp';
hNC    = hNC';
NCa    = NCa';
NCc = NCc';

%% check binomial distribution


NCsig = zeros(1,size(NCa,1));
NCbinom = zeros(1,size(NCa,1));
NCpBinom = zeros(1,size(NCa,1));

center = ceil(nBins/2);
checkWin = 500/statsBinWidth;

sigWin = 3; 

for ch = 1:size(NCa,1)
    
    h = hNC(ch,:);
    checkSig = mask2bounds(h);
    checkSig = abs(checkSig(:,2) - checkSig(:,1)) + 1;
    if max(checkSig) >= sigWin
        NCsig(ch) = 1;
    end
    
    s = sum(NCc(ch,center+1:center+checkWin)); n = sum(NCc(ch,[center-checkWin:center-1,center+1:center+checkWin])); 
    pout = myBinomTest(s,n,0.5,'two');
    NCpBinom(ch) = pout;
    NCbinom(ch) = s/n;
    
end

sig = NCsig;
pBinom = NCpBinom;
binom = NCbinom;

sig = logical(sig);

fprintf(['\nSig. Mod      : ',num2str(sum(sig)/length(sig)*100),' ',num2str(sum(sig)), ' / ',num2str(length(sig)), '\n'])
fprintf(['Sig. Sideness : ',num2str(sum(pBinom(sig) < 0.05) / sum(sig)*100), ' ', num2str(sum(pBinom(sig) < 0.05)), ' / ',num2str(sum(sig)), '\n'])


%% Plotting

colPath = './../data/';
load(fullfile(colPath,'col.mat'))

h1 = figure;

middle = movmean(-win+(plotBinWidth/2):plotBinWidth:win,2);
middle(1) = [];
for subj = 1:17
    subjii = subjID == subj;
    ii = NCsig & subjii;% 

    normFactor = NCn(:,:,logical(ii));
    sum(ii)
    normHCa = NCa(logical(ii),:) / mean(normFactor(:));
    data = smoothdata(mean(normHCa), 'movmean',1);
    b = plot(middle,data); hold on;
    b.LineWidth = 1.5;
    b.Color = col(subj,:);
    b.Color(4) = 0.4;  

end
fig = gcf;
fig.Color = 'w';
xlim([-1000 1000])
ax = gca;
ax.XTick = -1000:250:1000;
vline(0)
ylabel('NC-R centers [Prop of baseline]')
title(state)
fig = gcf;
fig.Color = 'w';


h2 = figure;

ii = NCsig;
normFactor = NCn(:,:,logical(ii));
normHCa = NCa(logical(ii),:) / mean(normFactor(:));
sem = std(normHCa) / sqrt(size(normHCa,1),'omitnan');
data = smoothdata(mean(normHCa), 'movmean',1);
[b, bPatch]= boundedline(middle,data,sem); hold on;
b.LineWidth = 2;
bPatch.FaceAlpha = 0.7;
null = squeeze(mean(normFactor,3)) / mean(normFactor(:));
sigLine  = plot(middle, quantile(null,0.99,1), '--'); hold on;
sigLine.Color = [0.5 0.5 0.5];     
sigLine.LineWidth = 2;
sigLine  = plot(middle, quantile(null,0.01,1), '--'); hold on;
sigLine.Color = [0.5 0.5 0.5];  
sigLine.LineWidth = 2;

fig = gcf;
fig.Color = 'w';

xlim([-1000 1000])
ax = gca;
ax.XTick = -1000:250:1000;

vline(0)

ylabel('NC-R centers [Prop of baseline]')
title(state)

fig = gcf;
fig.Color = 'w';
% savepdf(h,fullfile(exportDir,['PRTH_NCNC_',state,'_sigChan.pdf']));
% saveas(h,fullfile(exportDir,['PRTH_NCNC_',state,'_sigChan.fig']));
% savepdf(h,fullfile(exportDir,['PRTH_NCNC_',state,'_sigChan_subject.pdf']));
% saveas(h,fullfile(exportDir,['PRTH_NCNC_',state,'_sigChan_subject.fig']));
