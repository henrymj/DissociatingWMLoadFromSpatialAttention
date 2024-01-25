function jump = ppa(rawTS,rateAcq,thresh)
% Detect whether peak-to-peak amplitude exceeds thresh.
%
% Last modifed 01.29.2015 by Joshua Foster
% - inputs: rawTS, rateAcq, thresh

% dat = squeeze(erp.arf.trial.data(trial,:,:));
% rawTS = dat(chan,:);

winSize = 15; % size of test window, ms
winSize = round(winSize/rateAcq); % size of test window, samples
winStep = 50; % size of step between windows, ms
winStep = round(winStep/rateAcq); % size of step, samples

wInd = 1; jump = zeros(size(rawTS));

while 1
    % determine portion of rawTS to test
    wEnd = wInd + winSize; 
    window = wInd:wEnd;
    
    minPeak = min(rawTS(window)); % minimum peak amplitude
    maxPeak = max(rawTS(window)); % maximum peak amplitude
    
    p2p = abs(maxPeak-minPeak); % peak to peak amplitude
    
    if p2p > thresh
        jump(window) = 1;
    end
    
    wInd = wInd + winStep;
    
    if wInd + winSize > length(rawTS)
        break
    end
end