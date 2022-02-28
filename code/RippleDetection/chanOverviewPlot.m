
function hand = chanOverviewPlot(cleanedEvents, plotStats, ch, win, fs, subject, chan_labels)

analyticRB = nan(size(cleanedEvents.rippleband{ch}));
analyticHF = nan(size(cleanedEvents.rippleband{ch}));
analytic100200 = nan(size(cleanedEvents.rippleband{ch}));

for rip = 1:size(cleanedEvents.rippleband{ch},1)
    analyticRB(rip,:) = abs(hilbert(cleanedEvents.rippleband{ch}(rip,:)));
    analyticHF(rip,:) = abs(hilbert(cleanedEvents.highFreqBand{ch}(rip,:)));
    analytic100200(rip,:) = abs(hilbert(cleanedEvents.band100200{ch}(rip,:)));
    
end
meanLFP(ch,:) = mean(cleanedEvents.raw{ch});
semLFP(ch,:) =  std(cleanedEvents.raw{ch}) / sqrt(size(cleanedEvents.raw{ch},1)); % nansem
meanRB = mean(analyticRB);
meanHF = mean(analyticHF);
mean100200 = mean(analytic100200);

center = round(size(meanRB,2)/2);

hand = figure('Position', [316 30 1020 885], 'Visible', 'off');

subtightplot(10,4, 1:12, [0.02, 0.04])
% yyaxis right
% ln1 = plot(-win:win,meanRB, 'k-'); hold on;
% ln1.LineWidth = 1.5;
% ln2 = plot(-win:win,meanHF + ceil(max(meanRB)), 'r-'); hold on;
% ln2.LineWidth = 1.5;
% ln3 = plot(-win:win,mean100200, 'g-'); hold on;
% ln3.LineWidth = 1.5;
% ylabel('bandpass \muV')

% yyaxis left
boundedline(-win:win, meanLFP(ch,center-win:center+win), semLFP(ch,center-win:center+win)); hold on

vline(0)

% xlabel('time (ms)')
ylabel('LFP \muV')
xlim([-fs*1.5 fs*1.5])

density = length(cleanedEvents.goodRipples{ch})/sum(cleanedEvents.recordingLength) * fs * 60;
title(sprintf('%s, %s: LFP, %i events, %2.1f /min', ...
    subject, chan_labels{ch}, size(cleanedEvents.raw{ch},1), density))
ax = gca;
ax.YAxis.Color = 'k';
% legend([ln1, ln2, ln3], {'rippleband','200+Hz','100-200Hz'})

ax.XTickLabel = ax.XTick * (1/fs) * 1000; %xaxis in ms 
box on

subtightplot(10,4, 13:24, [0.1, 0.04])
winZoom = 0.25 * fs;

yyaxis right
ln1 = plot(-winZoom:winZoom,meanRB(center-winZoom:center+winZoom), 'k-'); hold on;
ln1.LineWidth = 2;
ln2 = plot(-winZoom:winZoom,meanHF(center-winZoom:center+winZoom) + ceil(max(meanRB)), 'r-'); hold on;
ln2.LineWidth = 2;
ln3 = plot(-winZoom:winZoom,mean100200(center-winZoom:center+winZoom), 'g-'); hold on;
ln3.LineWidth = 2;
ylabel('bandpass ampitude \muV')

yyaxis left
boundedline(-winZoom:winZoom, meanLFP(ch,center-winZoom:center+winZoom), semLFP(ch,center-winZoom:center+winZoom)); hold on

vline(0)

% xlabel('time (ms)')
ylabel('LFP \muV')
xlim([-winZoom winZoom])

ax = gca;
ax.YAxis(2).Color = 'k';
legend([ln1, ln2, ln3], {'rippleband','200+Hz','100-200Hz'})

ax.XTickLabel = ax.XTick * (1/fs) * 1000; %xaxis in ms 
xlabel('time to ripple center [ms]')
            
% 
subtightplot(10,4, [25,26,29,30], [0.04, 0.06])
hgram = histogram(plotStats.rippleAmp{ch}, 100, 'EdgeColor','None'); hold on;
meanVal = mean(plotStats.rippleAmp{ch});
stdVal = std(plotStats.rippleAmp{ch});
% xlabel(sprintf('Ripple Amp: %.2f (%.2f) [\\muV]', meanVal, stdVal))
xlabel('Ripple Amp [\\muV]')
plot([meanVal meanVal], [0 1000], 'k-');
plot([meanVal+stdVal meanVal+stdVal], [0 1000], 'r--');
plot([meanVal-stdVal meanVal-stdVal], [0 1000], 'r--');
ylim([0 max(hgram.Values)])
ylabel('ripple counts')

subtightplot(10,4, [27,28,31,32], [0.04, 0.06])
hgram = histogram(plotStats.duration{ch}, 'BinWidth', 1000/fs, 'EdgeColor','None'); hold on;
meanVal = mean(plotStats.duration{ch});
stdVal = std(plotStats.duration{ch});
% xlabel(sprintf('Ripple Dur: %.2f (%.2f) [ms]', meanVal, stdVal))
xlabel('Ripple Dur [ms]')
plot([meanVal meanVal], [0 1000], 'k-');
plot([meanVal+stdVal meanVal+stdVal], [0 1000], 'r--');
plot([meanVal-stdVal meanVal-stdVal], [0 1000], 'r--');
ylim([0 max(hgram.Values)])
ylabel('ripple counts')

subtightplot(10,4, [33,34,37,38], [0.04, 0.06])
hgram = histogram(plotStats.oscFreq{ch}, 100, 'EdgeColor','None'); hold on;
meanVal = mean(plotStats.oscFreq{ch});
stdVal = std(plotStats.oscFreq{ch});
xlabel('Ripple Osc Freq [Hz]')
plot([meanVal meanVal], [0 1000], 'k-');
plot([meanVal+stdVal meanVal+stdVal], [0 1000], 'r--');
plot([meanVal-stdVal meanVal-stdVal], [0 1000], 'r--');
ylim([0 max(hgram.Values)])
ylabel('ripple counts')

subtightplot(10,4, [35,36,39,40], [0.04, 0.06])

hgram = histogram(log10(plotStats.InterRipPeriod{ch}), 100, 'EdgeColor','None'); hold on;
ax = gca;
ticks = floor(min(ax.XTick)):ceil(max(ax.XTick));
ax.XTick = ticks;
ax.XTickLabels = 10.^ax.XTick;

meanVal = mean(plotStats.InterRipPeriod{ch});
stdVal = std(plotStats.InterRipPeriod{ch});
xlabel('Inter-Ripple Interval [secs]')
plot(log10([meanVal meanVal]), [0 1000], 'k-');
% plot(log10([meanVal+stdVal meanVal+stdVal]), [0 1000], 'r--');
% plot(log10([meanVal-stdVal meanVal-stdVal]), [0 1000], 'r--');
ylim([0 max(hgram.Values)])
ylabel('IRI counts')

fig = gcf;
fig.Color = [1 1 1];


return


