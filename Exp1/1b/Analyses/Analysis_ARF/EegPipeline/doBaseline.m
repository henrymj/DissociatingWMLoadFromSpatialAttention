function baselined = doBaseline(data,times,baselineStart,baselineEnd)
% Baseline segmented EEG
%
% INPUTS:
% data: matrix of segmented data to baseline (trials x channels x electrode)
% times: the sample timepoints for the segmented data (e.g., -500, -498, -496...., 1496, 1498, 1450)
% baselineStart: start of baseline in ms (e.g., -500 ms)
% baselineEnd: end of baseline period in ms (e.g., 0 ms)
% 
% OUTPUT: 
% baselined: the baselined data.
%
% Written by Joshua J. Foster, January 28, 2016

fprintf('baselining data... \n')

% get dimensions from data
nTrials = size(data,1); nChans = size(data,2); nSamps = size(data,3);

% create baseline index
bInd = ismember(times,baselineStart:baselineEnd); 

% preallocate array for baselined data
baselined = nan(size(data));

for t = 1:nTrials
    for chan = 1:nChans;
        dat = squeeze(data(t,chan,:));    % grab the raw time-series
        base = mean(dat(bInd));           % calcuate mean during baseline period
        baselined(t,chan,:) = dat-base;   % do baseline subtraction
    end
end

