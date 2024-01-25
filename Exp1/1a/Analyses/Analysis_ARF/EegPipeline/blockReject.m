function blocking = blockReject(rawTS,rateAcq)
% blockReject implements X-within-Y-of-peak rejection procedure as described
% in "An Introduction to the Event Related Potential Technique" by Steven J. Luck (pp. 168).

% Last modified 09.01.2015 by Kirsten Adam
% - Changed # of points (x) to be determined by sampling rate

% Modified 01.29.2016 by Joshua Foster
% - inputs are now acqRate and raw time-series

% give the relevent segment as input!
% dat = squeeze(erp.arf.trial.data(trial,:,:));

winSize = 200; % size of test window, ms
winSize = round(winSize/rateAcq); % size of test window in sample points.
winStep = 50; % size of step between windows, ms
winStep = round(winStep/rateAcq); % size of step between windows in sample points.

y = 1; % threshold - how close a value must be to the peak to be counted, microvolts.
x =  60 ./rateAcq; % # of points that must within Y of peak


wInd = 1; blocking = zeros(size(rawTS));

while 1
    % determine portion of raw TS to test.
    wEnd = wInd +winSize;
    window = wInd:wEnd;
    
    % calculate the min and max for the window
    maxPeak = max(rawTS(window));
    minPeak = min(rawTS(window));
    
    diffFromMax = maxPeak - rawTS(window); % always +ve
    diffFromMin = rawTS(window) - minPeak; % always +ve
    
    xwithinMax = diffFromMax < y;
    xwithinMin = diffFromMin < y;
    
    % mark the window as blocking if more than x points y within the min or
    % max
    if sum(xwithinMax) > x || sum(xwithinMin) > x
        blocking(window) = 1; % changed from 1 to chan so that we can visualize which chans were counted as blocking! 
    end

    wInd = wInd + winStep;
    
    if wInd + winSize > length(rawTS)
        break
    end   
end