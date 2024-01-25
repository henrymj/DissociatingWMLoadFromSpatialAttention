function eyeData = segment_eye(eyeData,preTime,postTime)

% Determine onset of each trial
eyeData.preTimeLockTime = preTime./eyeData.rateAcq;
eyeData.postTimeLockTime = postTime./eyeData.rateAcq;
eyeData.totalTime = eyeData.preTimeLockTime + eyeData.postTimeLockTime + 1;
% number time points
maxTrial = eyeData.preTimeLockTime+eyeData.postTimeLockTime+1;

% preallocate some matrices

data = ones(eyeData.nTrials,eyeData.nChans,maxTrial).*NaN; % Data matrix: #Conditions x #Trials/condition x length of trial x #Channels
tTimes = ones(eyeData.nTrials).*NaN;
tNums = ones(eyeData.nTrials).*NaN;
rejSaccMat = nan(1,eyeData.nTrials);
rejBlMat = nan(1,eyeData.nTrials);
% for i = 1:eyeData.cond.nConds
%fprintf('%d \t Percent Complete \n', round((i/eyeData.cond.nConds)*100))
dMat = ones(eyeData.nTrials,eyeData.nChans,maxTrial).*NaN;
% rejects = zeros(3,eyeData.cond.nTrials);
%     cc = eyeData.cond.condCodes(:,i);
%     cInd = ismember(eyeData.eventCodes,cc);
cInd = eyeData.startTimes;

cnt = 1; tCnt = 0; rCnt = 0;
for ii  = 1:length(eyeData.startTimes)
    %         if cnt >= nCodes % start counting after number of event codes per trial
    %             cCnt = sum(cInd(ii-nCodes+1:ii));
    %             if cCnt == nCodes
    tCnt = tCnt+1;
    % Store event times for each trial
    %                 tTimes(tCnt,:) = eyeData.sampleTimes(ii-nCodes+1:ii);
    %                 tNums(i,tCnt) = eyeData.trial.trialNums(ii-nCodes+1);
    %                 tCodes(i,tCnt) = eyeData.eventCodes(ii-nCodes+1);
    % Determine start and stop of trial
    tStart = eyeData.startTimes(ii)-eyeData.preTime; % in ms, subtract actual milleseconds!!!
    tEnd = eyeData.startTimes(ii)+eyeData.postTime; % in ms, subract actual milliseconds!
    %                 tLength = maxTrial; %tEnd-tStart+1;
    tWindow = tStart:tEnd;
    tWindowInd = ismember(eyeData.sampleTimes,tWindow); % index which times with a logical (will drop to actual number of samples after dividing by rate acq)
    
    if sum(tWindow<=0) ~= 0
        %%%%% if there is a negative number (not all times got
        %%%%% saved b/c switched on late--- make it a NaN!!! KA
        %%%%% 1/2015
        rawTS = NaN(eyeData.nChans,length(tWindow));
        
    else
        %%%% otherwise proceed as normal. KA 1/2015
        %                     rawTS = eyeData.data(eyeIndex,tWindow);
        %changed by Sutterer 1/25/16 this way if we for some reason switch
        %eyes during a session we won't grab the wrong eye. 
        leftEyeIndex = [1;0;1;0;1;0;1];%%%% always get the last one (distance)
        leftEyeIndex = logical(leftEyeIndex);
        rightEyeIndex = [0;1;0;1;0;1;1]; %%% % always get the last one (distance)
        rightEyeIndex = logical(rightEyeIndex);
        if eyeData.RecordedEye(ii) == 1
            rawTS = eyeData.all_data_gaze(leftEyeIndex,tWindowInd);
        else
            rawTS = eyeData.all_data_gaze(rightEyeIndex,tWindowInd);
        end
        %[eyeData.hx(tWindowInd); eyeData.hy(tWindowInd); eyeData.pa(tWindowInd)];
        % change the saved labels
        tLength = size(rawTS,2);
        if sum(ismember(tWindow,eyeData.blinkStartInd))>0|sum(ismember(tWindow,eyeData.blinkEndInd))>0
            rejBlMat(ii) = 1;
        else
            rejBlMat(ii) = 0;
        end
        
        if sum(ismember(tWindow,eyeData.saccStartInd))>0|sum(ismember(tWindow,eyeData.saccEndInd))>0
            rejSaccMat(ii) = 1;
        else
            rejSaccMat(ii) = 0;
        end
    end
    % Store time series data if raw data passes arf
    dMat(tCnt,1:eyeData.nChans,1:tLength) = rawTS;
    %                 cnt = 0; % restart counter
    %             end
    %         end
    cnt = cnt+1;
end
% Collect raw data and provide rejection stats
data(1:tCnt,1:eyeData.nChans,:) = dMat(1:tCnt,:,:);

eyeData.epoched = squeeze(data);
eyeData.parserRejBlink = rejBlMat; %logical of if there were blinks
eyeData.parserRejSacc = rejSaccMat;

end
% end