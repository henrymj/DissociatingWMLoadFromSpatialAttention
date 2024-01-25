function eMove = step_eyeTrack(xDat,yDat,winStep,winSize,thresh,rateAcq,viewDist)
% function detects eye movements

% xDat and yDat are vectors of values for a single trial

%thresh = erp.arf.thresh.eMove;
%winSize = 100; % size of test window, ms
%winStep = 50; % size of step between windows, ms
winSize = round(winSize/rateAcq); % size of test window, samples
winStep = round(winStep/rateAcq); % size of step, samples

rawX = xDat;
rawY = yDat;

wInd = 1; eMove = zeros(1,size(xDat,2));
while 1
    % determine portion of rawTS to test
    wEnd = wInd + winSize; 
    p = round(mean([wInd wEnd]));
    window = wInd:wEnd;
    
    prePeakX = mean(rawX(wInd:p)); % mean amplitude prior to pointer index
    postPeakX = mean(rawX(p:wEnd)); % mean amplitude after pointer index
       
    prePeakY = mean(rawY(wInd:p)); % mean amplitude prior to pointer index
    postPeakY = mean(rawY(p:wEnd)); % mean amplitude after pointer index
    
%     stepDeg = abs(postPeak-prePeak); % peak to peak amplitude
    stepDeg = pix2deg_changeInPos(prePeakX,postPeakX,prePeakY,postPeakY,viewDist);
    
    if stepDeg > thresh
        eMove(window) = 1; 
    end
    
    wInd = wInd + winStep;
    
    if wInd + winSize > length(rawX)
        break
    end
end