function eMove = step_t(rawTS,rateAcq,winStep,winSize,thresh)
% function detects eye movements
% edited by KA to accept single trial data rather than the full dataset ---
% save massive amounts of time!!!!! 

% dat = squeeze(erp.arfDat.data(trial,:,:));
% dat = rawTS;
%thresh = erp.arf.thresh.eMove;
%winSize = 100; % size of test window, ms
%winStep = 50; % size of step between windows, ms
winSize = round(winSize/rateAcq); % size of test window, samples
winStep = round(winStep/rateAcq); % size of step, samples

% rawTS = dat(chan,:);
wInd = 1; eMove = zeros(1,size(rawTS,2));
while 1
    % determine portion of rawTS to test
    wEnd = wInd + winSize; 
    p = round(mean([wInd wEnd]));
    window = wInd:wEnd;
    
    prePeak = mean(rawTS(wInd:p)); % mean amplitude prior to pointer index
    postPeak = mean(rawTS(p:wEnd)); % mean amplitude after pointer index
    
    stepAmp = abs(postPeak-prePeak); % peak to peak amplitude
    
    if stepAmp > thresh
        eMove(window) = 1; 
    end
    
    wInd = wInd + winStep;
    
    if wInd + winSize > length(rawTS)
        break
    end
end