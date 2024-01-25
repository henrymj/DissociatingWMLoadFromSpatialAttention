function erp = Rereference(erp)
% Rereference to the average of the Left and Right mastoids
% Created by David W. Sutterer 12/4/2015
% Checked by Kirsten C.S. Adam on: 12/14/2015
% Checked by Joshua J. Foster on: 1/26/2016
% Based on Introduction to the Event-Related Potential Technique pg. 107-108

chanLabels = erp.chanLabels;          % grab channel labels from erp struct
notIncluded = erp.droppedElectrodes;  % grab the list of electrodes to drop

leftMastoid = ismember(chanLabels,'TP9')';               % TP9 is the left mastoid, TP10 is the right mastoid (and the online reference)
l = length(notIncluded); notIncluded{l+1} = 'TP9';       % add the offline reference (TP9) to the notIncluded vector
notIncludedInd = ismember(chanLabels,notIncluded)';   
chanLabels_reref = chanLabels(~notIncludedInd);     
skipReRef = {'HEOG','VEOG','StimTrak'};                  % eye channels and Stimtrak are last and have their own references so don't rereference them
Chan2ReRef = sum(~ismember(chanLabels_reref,skipReRef)); % drop TP9 from our reReferenceSet and don't include eye channels and stim
tmpreRef = erp.data(~notIncludedInd,:); % a    
mastoidValue = erp.data(leftMastoid,:) ./ 2; % r/2
for chan = 1:Chan2ReRef %NOTE This for loop is set up for HEOG,VEOG, and Stimtrak (which have their own face reference) as the last 3 channels. If you change the order make sure you restructure the script to skip rereferencing those channels!
    tmpreRef(chan,:) = tmpreRef(chan,:) - mastoidValue; % a - (r/2)
end

%save chan labels for moving forward
erp.data = tmpreRef;
erp.chanLabels = chanLabels_reref;
erp.nChans = sum(~notIncludedInd);
erp.nArfChans = sum(~notIncludedInd)-1; % don't count Stimtrak as a chan to be looked at for arfs
erp.chanCoordinates = erp.chanCoordinates(~notIncludedInd);   % JJF: update this....
