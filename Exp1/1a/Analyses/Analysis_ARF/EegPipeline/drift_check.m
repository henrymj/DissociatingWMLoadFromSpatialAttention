function drift = drift_check(rawTS,rateAcq,thresh)
% check for extreme drift!!!
% just take an absolute voltage difference for the first quarter versus the
% last quarter of the time window of interest
% 
% Written by Kirsten Adam; July 11, 2016

drift = zeros(size(rawTS));

nSamples = length(rawTS);
breaks = round(linspace(1,nSamples,5)); % split epoch into quarters

firstQuarter = mean(rawTS(breaks(1):breaks(2)));
lastQuarter = mean(rawTS(breaks(4):breaks(5)));

if abs(firstQuarter-lastQuarter)>thresh
    drift = ones(size(rawTS));
end
