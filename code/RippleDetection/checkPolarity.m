% load ch x samples "raw"
function polCheck = checkPolarity(raw, sfreq) 

raw(isnan(raw)) = 0;
% detect positive and negative deflections
[SOmat, ~] = detectSO(raw, sfreq, [], 20);

% high gamma bandpass and find envelope
if sfreq/2 > 190
    [b,a] = butter(3, [70 190]/(sfreq/2));
else
    [b,a] = butter(3, 70/(sfreq/2), 'high');
end
    
HG_env = abs(hilbert(filtfilt(b,a,raw')))';

%% determine polarity by comparing the HG +/-win around positive vs. negative peaks
% if positive peaks have larger HG then polarity = 1

win = 100;
polCheck = NaN(size(raw,1),1);

for ch = 1:size(raw,1)
    
    pos = find(SOmat(:,3)==ch & SOmat(:,2)>0);
    neg = find(SOmat(:,3)==ch & SOmat(:,2)<0);
    
    pos(SOmat(pos,1)<=win | SOmat(pos,1)>size(raw,2)-win) = [];
    neg(SOmat(neg,1)<=win | SOmat(neg,1)>size(raw,2)-win) = [];

    p_HG = zeros(size(pos));
    n_HG = zeros(size(neg));
        
    for p = 1:numel(pos)
        p_HG(p) = mean(HG_env(ch,SOmat(pos(p),1)-win:SOmat(pos(p),1)+win));
    end
    
    for n = 1:numel(neg)
        n_HG(n) = mean(HG_env(ch,SOmat(neg(n),1)-win:SOmat(neg(n),1)+win));
    end
    
    if mean(p_HG) > mean(n_HG)
        polCheck(ch) = 1;
    else
        polCheck(ch) = -1;
    end
    
    fprintf(['channel ', num2str(ch), ' polarity : ', num2str(polCheck(ch)), '\n'])

end

% polCheck = logical(polCheck);

return
