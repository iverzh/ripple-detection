% AnalyzeCrossCorrelogram.m 
% Script used to produce all neocortical-hippocampal ripple cross correlogram
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


addpath(genpath('.'))


%% Inputs

subj_list_full = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
                  'S11','S12','S13','S14','S15','S16','S17'};

mode           = 'HC'; %HC, SSR or SWR
state          = 'wake'; %wake or sleep
runList        = 1:17; %list of subjects you want to run

PRTHfolderHC   =  ['./../data/CrossCorrelogram/',state,'/',mode]; %Preprocessed crosscorrelogram data
badChanDir     = './../data/bad_chan';
exportDir      =     './../Figures';
if ~isfolder(exportDir); mkdir(exportDir); end

FDRwin        = 500; % +- window (in ms) to compute FDR significance
statsBinWidth = 25;  %histogram bin width (in ms) for computing all stats (significance, sidedness)
plotBinWidth  = 25; %histogram bin width (in ms) for plotting
flip          = false; %flip cross correlogram axes
fs            = 1000; %sampling rate

%% load peri-ripple time histogram data

HCp = [];
HCpFDR = [];
HCa = [];
HCc = [];
HCn = [] ; 

HCa_plot = [];
HCn_plot = [];

HCa_parcel = [];
HCa_parcel_plot = [];

subjNames = {};
subjID = [];
parcelLocs = [];
c = 0;
nRip = 0;
for subj = runList 
    subject = subj_list_full{subj};
    
    if flip
        HC = load(fullfile(PRTHfolderHC,[subject,'_',mode,'_PRTH_flip.mat']));
    else
        HC = load(fullfile(PRTHfolderHC,[subject,'_',mode,'_PRTH.mat']));
    end
    
    % compute histogram parameters
    nIter = size(HC.subjPRTH.(mode).nullPRTH{1,1},1);
    win = ceil(size(HC.subjPRTH.(mode).nullPRTH{1,1},2)/2);
    edges = -win+(statsBinWidth/2):statsBinWidth:win;
    edgesPlot = -win+(plotBinWidth/2):plotBinWidth:win;
    times = -win:win;
    nBins = 2*win/statsBinWidth - 1;
    nBinsPlot = 2*win/plotBinWidth - 1;
    center = ceil(nBins/2);

    load(fullfile(badChanDir,[subject,'_bad_chan.mat']))    
    
    for ch = 1:length(HC.subjPRTH.chanLabels)
        chanLabel = HC.subjPRTH.chanLabels{ch};
        
        if ~contains(chanLabel,bad_chan) 
            for HCch = 1:length(HC.subjPRTH.chanLabelsHC)
                HCchanLabel = HC.subjPRTH.chanLabelsHC{HCch};
                if ~contains(HCchanLabel,bad_chan)
                    
                    countsTemp = smoothdata(HC.subjPRTH.(mode).eventPRTH{ch,HCch}, 'gaussian',250);
                    nullTemp   = smoothdata(HC.subjPRTH.(mode).nullPRTH{ch,HCch}, 2, 'gaussian',250);
                    countsOrig = zeros(1,nBins);
                    nullOrig = zeros(nIter,nBins);

                    for e = 1:length(edges)-1
                       ii = times > edges(e) & times < edges(e+1);
                       countsOrig(e) = sum(countsTemp(ii));
                       nullOrig(:,e) = sum(nullTemp(1:nIter,ii),2);
                    end

                    counts = countsOrig; 
                    null   = nullOrig; 

                    countsTemp = smoothdata(HC.subjPRTH.(mode).eventPRTH{ch,HCch}, 'gaussian',50);
                    nullTemp   = smoothdata(HC.subjPRTH.(mode).nullPRTH{ch,HCch}, 2, 'gaussian',50);
                    countsOrig = zeros(1,nBinsPlot);
                    nullOrig = zeros(nIter,nBinsPlot);

                    for e = 1:length(edgesPlot)-1
                       ii = times > edgesPlot(e) & times < edgesPlot(e+1);
                       countsOrig(e) = sum(countsTemp(ii));
                       nullOrig(:,e) = sum(nullTemp(1:nIter,ii),2);
                    end

                    countsPlot = countsOrig; %smoothdata(countsOrig, 'gaussian',5);
                    nullPlot   = nullOrig; %smoothdata(nullOrig, 2, 'gaussian',5);


                    n =  HC.subjPRTH.(mode).nRipples{ch,HCch};
                    switch mode
                        case {'HC','SSR','SWR'}
                            HCp = [HCp, findPercentile(null,counts,'pval')];
                            HCpFDR = [HCpFDR, findPercentile(null(:,center-(FDRwin/statsBinWidth):center+(FDRwin/statsBinWidth)),counts(center-(FDRwin/statsBinWidth):center+(FDRwin/statsBinWidth)),'pval')];
                        case 'none'
                            HCp = [HCp, findPercentile(null,counts,'pval2side')];  
                            HCpFDR = [HCpFDR, findPercentile(null(:,center-(FDRwin/statsBinWidth):center+(FDRwin/statsBinWidth)),counts(center-(FDRwin/statsBinWidth):center+(FDRwin/statsBinWidth)),'pval2side')];

                    end


                    HCa = [HCa, counts/(n*statsBinWidth)*fs]; 
                    HCa_plot = [HCa_plot, countsPlot/(n*plotBinWidth)*fs]; 

                    HCc = [HCc, HC.subjPRTH.(mode).eventPRTH{ch,HCch}]; 

                    HCn = [HCn, null/(n*statsBinWidth)*fs]; 
                    HCn_plot = [HCn_plot, nullPlot/(n*plotBinWidth)*fs]; 


                    c = c + 1;
                    subjNames{c} = [subject, chanLabel];
                    subjID(c) = subj;

                end
                
                
            end
        end
        
    end
    
    fprintf([subject, '\n'])
    
end


%% Compute FDR

FDRbins = 2*FDRwin/statsBinWidth+1;

[hHC, ~, ~, adjHCp]  = fdr_bh(HCpFDR, 0.05, 'pdep','yes');

adjHCp = reshape(adjHCp, [length(adjHCp)/FDRbins, FDRbins]);
hHC    = reshape(hHC,    [FDRbins, length(hHC)/FDRbins]);

HCa = reshape(HCa, [nBins, length(HCa)/nBins]);
HCa_plot = reshape(HCa_plot, [nBinsPlot, length(HCa_plot)/nBinsPlot]);

HCc = reshape(HCc, [2*win-1, length(HCc)/(2*win-1)]);

HCn = reshape(HCn, [nIter, nBins, length(HCn)/nBins]);
HCn_plot = reshape(HCn_plot, [nIter, nBinsPlot, length(HCn_plot)/nBinsPlot]);

adjHCp = adjHCp';
hHC    = hHC';
% 
HCa    = HCa';
HCa_plot    = HCa_plot';

HCc = HCc';

%% check binomial distribution

HCsig = zeros(1,size(HCa,1));
HCbinom = zeros(1,size(HCa,1));
HCpBinom = zeros(1,size(HCa,1));

checkWin = 500;
center = ceil((2*win+1)/2);
sigWin = 3; 

for ch = 1:size(HCa,1)

    h = hHC(ch,:);
    checkSig = mask2bounds(h); 
    checkSig = abs(checkSig(:,2) - checkSig(:,1)) + 1;
    if max(checkSig) >= sigWin
        HCsig(ch) = 1;
    end
    
    s = sum(HCc(ch,center+1:center+checkWin)); n = sum(HCc(ch,[center-checkWin:center-1,center+1:center+checkWin])); 
    pout = myBinomTest(s,n,0.5,'two');
 
    HCpBinom(ch) = pout;
    HCbinom(ch) = s/n;
    
end


sig = HCsig;
pBinom = HCpBinom;
binom = HCbinom;

sig = logical(sig);

fprintf(['Sig. Mod      : ',num2str(sum(sig)/length(sig)*100),' ',num2str(sum(sig)), ' / ',num2str(length(sig)), '\n'])
fprintf(['Sig. Sideness : ',num2str(sum(pBinom(sig) < 0.05) / sum(sig)*100), ' ', num2str(sum(pBinom(sig) < 0.05)), ' / ',num2str(sum(sig)), '\n'])

if flip
    fprintf(['NC-R Leading  : ',num2str(sum(binom(sig & pBinom < 0.05) < 0.5) /sum(pBinom(sig) < 0.05)*100), ' ', num2str(sum(binom(sig & pBinom < 0.05) < 0.5)), ' / ',num2str(sum(pBinom(sig) < 0.05)), '\n'])
else
    fprintf(['NC-R Leading  : ',num2str(sum(binom(sig & pBinom < 0.05) > 0.5) /sum(pBinom(sig) < 0.05)*100), ' ', num2str(sum(binom(sig & pBinom < 0.05) > 0.5)), ' / ',num2str(sum(pBinom(sig) < 0.05)), '\n'])
end
%%
colPath = './../data/';
load(fullfile(colPath,'col.mat'))

close all
h1 = figure;

middle = movmean(-win+(plotBinWidth/2):plotBinWidth:win,2);
middle(1) = [];
for subj = 1:17
    subjii = subjID == subj;
    ii = HCsig & subjii;

    normFactor = HCn_plot(:,:,logical(ii));
    normHCa = HCa_plot(logical(ii),:) / mean(normFactor(:));
    data = smoothdata(mean(normHCa), 'movmean',1);
    b = plot(middle,data); hold on;
    b.LineWidth = 1.5;
    b.Color = col(subj,:);
    b.Color(4) = 0.4;


end

xlim([-750 750])
ax = gca;
ax.XTick = -750:250:750;
% switch state
%         case 'sleep'
%             ylim([0.8 1.6])
% end

vline(0)
ylabel([mode,' centers [Prop of baseline]'])
xlabel('NC-R centers [ms]')
title(state)

fig = gcf;
fig.Color = 'w';

h2 = figure;
ii = HCsig;
normFactor = HCn_plot(:,:,logical(ii));
    
normHCa = HCa_plot(logical(ii),:) / mean(normFactor(:));
sem = std(normHCa, 'omitnan') / sqrt(size(normHCa,1));
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

xlim([-750 750])
ax = gca;
ax.XTick = -750:250:750;
% switch state
%         case 'sleep'
%             ylim([0.8 1.6])
% end

vline(0)
ylabel([mode,' centers [Prop of baseline]'])
xlabel('NC-R centers [ms]')
title(state)

fig = gcf;
fig.Color = 'w';

% if flip
%     savepdf(h1,fullfile(exportDir,['PRTH_',mode,'_',state,'_binWidth_',num2str(plotBinWidth),'_sigChan_flip_subject.pdf']))
%     saveas(h1,fullfile(exportDir,['PRTH_',mode,'_',state,'_binWidth_',num2str(plotBinWidth),'_sigChan_flip_subject.fig']))%
%     savepdf(h2,fullfile(exportDir,['PRTH_',mode,'_',state,'_binWidth_',num2str(plotBinWidth),'_sigChan_flip.pdf']))
%     saveas(h2,fullfile(exportDir,['PRTH_',mode,'_',state,'_binWidth_',num2str(plotBinWidth),'_sigChan_flip.fig']))
% else
%     savepdf(h1,fullfile(exportDir,['PRTH_',mode,'_',state,'_binWidth_',num2str(plotBinWidth),'_sigChan_subject.pdf']))
%     saveas(h1,fullfile(exportDir,['PRTH_',mode,'_',state,'_binWidth_',num2str(plotBinWidth),'_sigChan_subject.fig']))
%     savepdf(h2,fullfile(exportDir,['PRTH_',mode,'_',state,'_binWidth_',num2str(plotBinWidth),'_sigChan.pdf']))
%     saveas(h2,fullfile(exportDir,['PRTH_',mode,'_',state,'_binWidth_',num2str(plotBinWidth),'_sigChan.fig']))
% end










