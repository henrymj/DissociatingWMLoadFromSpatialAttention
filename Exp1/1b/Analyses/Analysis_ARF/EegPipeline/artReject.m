function erp = artReject(erp)
% function for rejecting trials
%
% Last modified 01.29.2016 by Joshua Foster
% - stripped out all mention of conds and removed blink and eye movement
% detection'
%
% Modified 2/25/16 : put eye blink detection back in because it's pretty reliable and super useful!! so much less clicking. KA  
% Modified 8/29/16 : edited drift checks, rewrote step function as step_t
% to make it a lot faster; added step function to check main cap channels
% for sudden dropout

fprintf('rejecting artifacts... \n')

nTrials = erp.arfDat.nTrials;
chanLabels = erp.chanLabels;
nChans = erp.nChans;
chanInd = 1:nChans;

% Preallocate matrices 
erp.arf.blink = NaN(1,nTrials);
erp.arf.eMove = NaN(1,nTrials);
erp.arf.artifactInd = NaN(1,nTrials);
erp.arf.grand = zeros(nChans,nTrials);      % matrix for labeling alllll the artifacts! 

%%%% temporary holders so that we can do a parfor loop!!! 
rateAcq = erp.rateAcq;
mark_blockingFull = NaN(nChans,nTrials);
mark_noiseFull = NaN(nChans,nTrials);
mark_driftFull = NaN(nChans,nTrials);
mark_dropoutFull = NaN(nChans,nTrials);

%%%% temporary holders for all  threshodl names!!! otherwise, the whole erp
%%%% structure will get shuttled over to the arf functions and this results
%%%% in hellllla slow behavior. KA 
noiseThr = erp.arf.noiseThr;
driftThr = erp.arf.driftThr;
dropoutWin = erp.arf.dropoutWin;
dropoutStep = erp.arf.dropoutStep;
dropoutThr = erp.arf.dropoutThr;
blinkWin = erp.arf.blinkWin;
blinkStep = erp.arf.blinkStep; 
blinkThr = erp.arf.blinkThr; 
eMoveWin = erp.arf.eMoveWin; 
eMoveStep = erp.arf.eMoveStep; 
eMoveThr = erp.arf.eMoveThr; 

fprintf('checking for blocking, noise, and drift... \n')


tic
for i = 1:nChans
    
    chanDat = squeeze(erp.arfDat.data(:,i,:));
    
    parfor t = 1:nTrials

        % get raw time-series for the trial and electrode
%         rawTS = squeeze(erp.arfDat.data(t,i,:));
        rawTS = chanDat(t,:);
        
       
        % check for blocking in all channels except StimTrak (and VEOG /
        % HEOG)
        checkChannel = ~ismember(chanLabels,'StimTrak');  % specify names of channels to skip
        if checkChannel(i)        
            block = blockReject(rawTS,rateAcq);
            mark = artDetect(block);
            mark_blockingFull(i,t) = mark;
        end
        
        % check for noise in all scalp channels
        checkChannel = ~ismember(chanLabels,{'HEOG','VEOG','StimTrak'}); % specify names of channels to skip        
        if checkChannel(i)
            noise = ppa(rawTS,rateAcq,noiseThr);
            mark = artDetect(noise);
            mark_noiseFull(i,t) = mark;
        end
        
        % check for extreme drift in all scalp channles 
        checkChannel = ~ismember(chanLabels,{'HEOG','VEOG','StimTrak'}); % specify names of channels to skip        
        if checkChannel(i)
            drift = drift_check(rawTS,rateAcq,driftThr);
            mark = artDetect(drift);
            mark_driftFull(i,t) = mark;
        end
        
       
        % check for extreme channel drop out (step function) in all scalp channles
        checkChannel = ~ismember(chanLabels,{'HEOG','VEOG','StimTrak'}); % specify names of channels to skip
        if checkChannel(i)
            dropout = step_t(rawTS,rateAcq,dropoutStep,dropoutWin,dropoutThr);
            mark = artDetect(dropout);
            mark_dropoutFull(i,t) = mark; 
        end
        
    end
end
%%%% put everything back outside of the parfor loop
erp.arf.blockingFull = mark_blockingFull;
erp.arf.noiseFull = mark_noiseFull;
erp.arf.driftFull = mark_driftFull;
erp.arf.dropoutFull = mark_dropoutFull;
toc
fprintf('checking for blinks... \n')

%%%%% check for blinks 
tic
mark = NaN(1,nTrials);
veogDat = squeeze(erp.arfDat.data(:,ismember(chanLabels,'VEOG'),:));
parfor t = 1:nTrials
    rawTS = veogDat(t,:);
    % Check for blinks using peak-to-peak amplitude
    blink = step_t(rawTS,rateAcq,blinkStep,blinkWin,blinkThr); % 33 = VEOG channel
    mark(t) = artDetect(blink);
end
erp.arf.blink = mark;
toc
fprintf('checking for eye movements... \n')

tic

markH = NaN(1,nTrials);
heogDat = squeeze(erp.arfDat.data(:,ismember(chanLabels,'HEOG'),:));
%%%%% check for horizontal eye movements 
parfor t = 1:nTrials
        rawTS = heogDat(t,:);
        eMoveH = step_t(rawTS,rateAcq,eMoveStep,eMoveWin,eMoveThr); % 33 = VEOG channel        
        markH(t) = artDetect(eMoveH);  
end

mark = markH;
erp.arf.eMove = mark;
toc

% create a vector "noise" that marks each trials as noise or no noise.
noiseSum = nansum(erp.arf.noiseFull,1);
noiseSum = squeeze(noiseSum);
nn = 0;
erp.arf.noise = ~ismember(noiseSum,nn);

% create a vector "noise" that marks each trials as noise or no noise.
blockSum = nansum(erp.arf.blockingFull,1);
blockSum = squeeze(blockSum);
nn = 0;
erp.arf.blocking = ~ismember(blockSum,nn);

% create a vector "drift" that marks each trials as noise or no noise.
driftSum = nansum(erp.arf.driftFull,1);
driftSum = squeeze(driftSum);
nn = 0;
erp.arf.drift = ~ismember(driftSum,nn);

% create a vector "dropout' that marks each trial as having channel drop
% out or not! 
dropoutSum = nansum(erp.arf.dropoutFull,1);
dropoutSum = squeeze(dropoutSum);
nn = 0; 
erp.arf.dropout = ~ismember(dropoutSum,nn);

% loop through trials and create an index of all artifacts
for t = 1:nTrials
    erp.arf.artifactInd(t) = erp.arf.blocking(t) | erp.arf.noise(t) | erp.arf.blink(t) | erp.arf.eMove(t) | erp.arf.drift(t) | erp.arf.dropout(t);
end

% save rejection statistics
erp.arf.totalArtProp = (sum(sum(erp.arf.artifactInd))/(nTrials)).*100;
erp.arf.blockingProp = (sum(sum(erp.arf.blocking))/(nTrials)).*100;
erp.arf.noiseProp = (sum(sum(erp.arf.noise))/(nTrials)).*100;
erp.arf.blinkProp =  (sum(sum(erp.arf.blink))/(nTrials)).*100;
erp.arf.eMoveProp = (sum(sum(erp.arf.eMove))/(nTrials)).*100;
erp.arf.driftProp =  (sum(sum(erp.arf.drift))/(nTrials)).*100;
erp.arf.dropoutProp =  (sum(sum(erp.arf.dropout))/(nTrials)).*100;

% print proportion of trials lost due to each kind of artifact.
fprintf('%d \tPercent trials rejected (total) \n', (round(erp.arf.totalArtProp)));
fprintf('%d \tPercent blocking \n', round(erp.arf.blockingProp));
fprintf('%d \tPercent noise \n', round(erp.arf.noiseProp));
fprintf('%d \tPercent eye movements \n', round(erp.arf.eMoveProp));
fprintf('%d \tPercent blinks \n', round(erp.arf.blinkProp));
fprintf('%d \tPercent drift \n', round(erp.arf.driftProp));
fprintf('%d \tPercent chan dropout \n', round(erp.arf.dropoutProp));

end
