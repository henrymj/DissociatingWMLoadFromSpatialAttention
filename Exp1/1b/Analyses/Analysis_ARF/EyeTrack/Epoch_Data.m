function eyeData = Epoch_Data(sn,dRoot)
%% load file 
subjects = sn; 
fname = [dRoot,num2str(subjects),'_PGSB.mat'];
load(fname); 
% fieldname = 'eye';
%.(fieldname) syntax allows you to access fields from a structured array
% RecordedEye = eye.RECORDINGS.('eye')% 1 = left 2 = right
RecordedEyeVec = [eye.RECORDINGS(:).('eye')];% so we can check if it changed over the course of the experiment (glasses and whatnot)
startEye = [1:2:length(RecordedEyeVec)]; %lets us just index the start of the recording so match up eye recording to trial start

%%
% FEVENT = contains timepoints with events (messages, eye blinks etc)
% FSAMPLE = contains all the sample we  actually want ! 

% eye.FEVENT.message = list of strings with messages (i.e. TRIALSTART)
% eye.FEVENT.eye = whether or not the line of text has anything to do with 
% eye.FEVENT.sttime = timestamp of each event 
% set some things up! 
eyeData = struct(); 
eyeData.rateAcq = 1; % 500 hz = sample every 2 ms 
eyeData.preTime = 400;  % absolute value in ms!! 
eyeData.postTime = 1550; 

%%%% settings for doing our artifact rejection
eyeData.winSize = 80; % ms --- size of sliding window that is looking for saccades
eyeData.stepSize = 10;  % ms ---- how much the window moves by 
eyeData.maxDeg = .5; % degrees of visual angle! if it's bigger than this reject it 

%%% Tell it the mode!! 1 = desktop, 2 = remote mode
eyeData.mode = 1; 
%%%% IF chin-rest mode, indicate distnce from the EYE to the SCREEN (in cm)
eyeData.DistanceToScreen = 73.8; % Measured on 9/21/16 - Don't move eyetracker from tape!
%%%%% IF remote mode, indicate distance from CAMERA to SCREEN (in cm) (distance
%%%%% from eye to camera is recorded in the data file and we will use that.
eyeData.DistanceToCamera = 37.6; % Measured on 9/21/16 - Don't move eyetracker from tape!



eyeData.RecordedEye=RecordedEyeVec;

% concatenate things so they are sensible to index and work with!! (at
% least in a way i understand)
eyeData.messageInd = {eye.FEVENT(:).message}; % look for messages we actually sent related to trial
eyeData.moveInd = {eye.FEVENT(:).codestring}; % look for start and beginning of saccades etc ! 
    % STARTSACC / ENDSACC; STARTFIX / ENDFIX ; STARTBLINK / ENDBLINK! 
eyeData.eventTimes = [eye.FEVENT(:).sttime];

% Get list of times associated with trial start and eye movement starts! 
%stuff we have to modify for our experiment
startTimeInd = strcmp('StimuliArray',eyeData.messageInd);
blinkStartInd = strcmp('STARTBLINK ',eyeData.moveInd); %space here is necessary
blinkEndInd = strcmp('ENDBLINK',eyeData.moveInd);
saccStartInd = strcmp('STARTSACC',eyeData.moveInd);
saccEndInd = strcmp('ENDSACC',eyeData.moveInd);
% get the index of times associated with cue onset
eyeData.startTimes = eyeData.eventTimes(startTimeInd);
%Added by Dave 1/15 to include parser data
eyeData.blinkStartInd = eyeData.eventTimes(blinkStartInd);
eyeData.blinkEndInd = eyeData.eventTimes(blinkEndInd);
eyeData.saccStartInd = eyeData.eventTimes(saccStartInd);
eyeData.saccEndInd = eyeData.eventTimes(saccEndInd);

% Epoch the data into trials 
eyeData.sampleTimes = [eye.FSAMPLE(:).time];
    % gaze coordinates 
    
    %at present I haven't spent the time to work in dynamic updating (in
    %case we switched eyes mid run). Shouldn't take too long though...
eyeData.gx = [eye.FSAMPLE(:).gx]; % gaze referenced x coords %  hx = head, gx = gaze 
eyeData.gy = [eye.FSAMPLE(:).gy]; % gaze referenced y coords
    % head referenced coordinates 
eyeData.hx = [eye.FSAMPLE(:).hx]; % head referenced x coords %  hx = head, gx = gaze 
eyeData.hy = [eye.FSAMPLE(:).hy]; % head referenced y coords
    % Pupil size / area 
eyeData.pa = [eye.FSAMPLE(:).pa]; % head referenced pupil size / area
  % distance!!!! (in mm, convert to cm!!)
eyeData.dist = (eye.FSAMPLE(:).hdata(3,:))./100; %%% 3rd row is distance; divide by 100 to scale to cm %%% only a meaningful variable in remote mode

% make an all data structure for the epoching routine!! 
eyeData.all_data_gaze= [eyeData.gx; eyeData.gy; eyeData.pa; eyeData.dist]; 
eyeData.all_data_head = [eyeData.hx; eyeData.hy; eyeData.pa; eyeData.dist]; 
eyeData.nTrials = sum(startTimeInd); %adds up logical index
%of ones and zeros where 1 indicates each trials start)

%%%% Dimension labels before just extracting 1 eye we recorded 
% 1:2 = left and right X, 3:4 = left and right Y, 5:6 = left and right pupil size, 7 =
% distance

%%%% Dimension labels after after:
eyeData.nChans = size(eyeData.all_data_gaze,1)-3; 


% add in a new field with epoched data
eyeData = segment_eye(eyeData,eyeData.preTime,eyeData.postTime); 

% do artifact rejection for this subject!!!!
[blinks,eMoves] = calc_eyeMove_dist_eyeTrack(eyeData);
eyeData.calculatedBadEyeVals = blinks;
eyeData.calculatedEMove = eMoves;

% print some information

fprintf(sprintf('Sub %d: Parser Blink Rate: %.2f \n',sn,sum(eyeData.parserRejBlink)./length(eyeData.parserRejBlink)))
fprintf(sprintf('Sub %d: Parser Saccade Rate: %.2f \n',sn,sum(eyeData.parserRejSacc)./length(eyeData.parserRejSacc)))
fprintf(sprintf('Sub %d: Calculated nonsensical position vals: %.2f \n',sn,sum(blinks)./length(blinks)))
fprintf(sprintf('Sub %d: Calculated eye move: %.2f \n',sn,sum(eMoves)./length(eMoves)))



