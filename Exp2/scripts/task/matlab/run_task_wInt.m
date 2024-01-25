clear;
close all;
Screen('CloseAll');
KbName('UnifyKeyNames');

%%  - EEG Session Vars
% screen - 1 for testing room
screenID = 0; %1;

% EEG porting
porting = false;
portCode = hex2dec('D050');

% eyetracking
eyes.tracking = false;
eyes.mode = 0; % 0 = chin rest, 1 = remote

%% - Setup & Task Details

% participant Info
prompt = {'SID', 'Participant Number', 'Participant Age', 'Participant Gender'};
defAns = {'999999', '999','99','X'};
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);
participant.ID = char(box(1));
participant.num = char(box(2));
participant.age = char(box(3));
participant.gender = char(box(4));

% Path details
sep = filesep;
paths.baseDir = fileparts(pwd);
paths.behavDir = [paths.baseDir, sep, 'DATA-Behavior'];
paths.eegDir = [paths.baseDir, sep, 'DATA-EEG', sep, participant.num];
if ~exist(paths.eegDir, 'dir')
   mkdir(paths.eegDir)
end
paths.fileBaseMain = ['p-' num2str(participant.num) '_task-changeDetection_Exp1'];  % use to get stimulus key mappings
paths.fileBase = ['p-' num2str(participant.num) '_task-changeDetection_Exp1_wInt'];

% Equipment Parameters
runPriority = 1;                   
equipment.viewDist_MM = 800;           % ! MEASURE !
equipment.pixPerMM = 3.6;              % ! MEASURE ! using MeasureDpi function

equipment.baseColor = [.5, .5, .5];
equipment.gammaVals = [1, 1, 1];

% Timing Parameters (in secs)               
timing.minITI = .55;
timing.maxITI = .95;
timing.baseline = .25;
timing.stepITI = .01;
timing.memoryArray = .25;
timing.postMemoryDelay1 = .5;
timing.postMemoryDelayInterupt = .15;
timing.postMemoryDelay2 = .85;
timing.postBreak = 1; 
timing.practiceFeedback = 1;
timing.ITIinterruptionWait = 1;

% Stimulus Details

stimulus.pInterupt = 1;

stimulus.fixationSize_DVA = .15;        % DVA = degrees of visual angle
stimulus.fixationColour = 0;
stimulus.regionHeight_DVA = 12;

stimulus.minEccentricity_DVA = 1.75;
stimulus.maxEccentricity_DVA = 5.25;
stimulus.minDotSize_DVA = .25;
stimulus.maxDotSize_DVA = .35;

stimulus.nBins = 8;
stimulus.nSlices = 40;
stimulus.nRows = 10;
stimulus.binSizes = [1, 3];
stimulus.minDots = 12;
stimulus.max1CloudDots = 48;
stimulus.max2CloudDots = 24;


stimulus.nTrialsInterupt = 192;
stimulus.nBlocks = 2;
stimulus.nTrialsPerBlock = stimulus.nTrialsInterupt/stimulus.nBlocks;

stimulus.sliceThickness = 2*pi / stimulus.nSlices; % radians per slice
stimulus.nDotsPerBin = stimulus.nSlices*stimulus.nRows/stimulus.nBins;
stimulus.binStarts = 1:stimulus.nDotsPerBin:stimulus.nDotsPerBin*(stimulus.nBins-1) + 1;

stimulus.colors = {'blue', 'green'};
stimulus.dotColors.blue =    [0, 0, 255]/255;
stimulus.dotColors.green =   [0, 60, 0]/255;
stimulus.dotColors.grey =   [35, 35, 35]/255;

stimulus.placeholderColor = equipment.baseColor;
stimulus.fixColor = 1;  

% Grab old key pairings if they exist; generate new ones if they don't.
oldInfoPath = [paths.fileBaseMain, '.mat'];
oldDir = cd(paths.behavDir);
if isfile(oldInfoPath)
    oldInfo = load(oldInfoPath);
    stimulus.keyCodes = oldInfo.stimulus.keyCodes; 
else
    stimulus.keyCodes = randsample([KbName('LeftArrow'), KbName('RightArrow')], 2);
end
cd(oldDir);

% Text parameters
text.instructionFont = 'Menlo';
text.instructionPoints = 16;        
text.instructionStyle = 0;          
text.instructionWrap = 80;
text.color = 1;

% Calculate Pixel Space 
equipment.MMperDeg = (equipment.viewDist_MM/2)*tan(deg2rad(2*stimulus.maxEccentricity_DVA))/stimulus.maxEccentricity_DVA;
equipment.PixPerDeg = equipment.pixPerMM*equipment.MMperDeg;

stimulus.fixationSize_pix = round(stimulus.fixationSize_DVA*equipment.PixPerDeg);
stimulus.regionHeight_pix = round(stimulus.regionHeight_DVA*equipment.PixPerDeg);

stimulus.minEccentricity_pix = round(stimulus.minEccentricity_DVA*equipment.PixPerDeg);
stimulus.maxEccentricity_pix = round(stimulus.maxEccentricity_DVA*equipment.PixPerDeg);

stimulus.minDotSize_pix = round(stimulus.minDotSize_DVA*equipment.PixPerDeg);
stimulus.maxDotSize_pix = round(stimulus.maxDotSize_DVA*equipment.PixPerDeg);

% Port Codes
if porting
    
    config_io;
    
    % trial codes
    %%% 1:96
    
    % block codes
    %%% 101:116 
    
    % Condition codes
    %%% 211:223
    %%% 200 + 
    %%% 30[= interruption trials] +
    %%% 1X = 1 cloud
    %%% 2X = 2 clouds
    %%% X1 = cloud(s) 1 bin width
    %%% X2 = clouds different widths
    %%% X3 = cloud(s) 3 bin width
    
    % Trial section codes
    
    %%% 150 (ITI starts)
    
    %%% trial code (Trial initiated - baseline starts)
    %%% condition code (Stim Display)
   
    %%% 160 (Delay Part 1 Starts)
    %%% 161 (Delay Interuption Starts - Interuption)
    %%% 162 (Delay Interuption Starts - no Interuption)
    %%% 163 (Delay Part 2 Starts)
    %%% 170 (Probe Display - no change)
    %%% 171 (Probe Display - change)
    %%% 180 (Response made)
    
    % BEGINNING AND END
    %%% 190 - BEGIN
    %%% 191 - END
    
end

%% - Trial Setup

% Build up Trial Details
rng('default');             
rng('shuffle');

% INTERRUPTION TRIALS
% Generated totally at random
interuptData.rng = rng;
interuptData.ITIs = datasample(timing.minITI:timing.stepITI:timing.maxITI, stimulus.nTrialsInterupt);

% Balance load
nClouds = [];  % balance within each block
for b = 1:stimulus.nBlocks
    nClouds = [nClouds; randsample([zeros(stimulus.nTrialsPerBlock/2,1) + 1; zeros(stimulus.nTrialsPerBlock/2,1) + 2], stimulus.nTrialsPerBlock)]; % half load 1, half load 2
end
interuptData.nClouds = nClouds;

for t = 1:stimulus.nTrialsInterupt
    interuptData.targetCloud(t).bins = get_bins_from_loc_width( ...
        randsample(1:8, 1), ... % loc
        randsample([1, 3], 1), ... % width
        stimulus.nBins); 
    
    target_color = randsample([1, 2], 1);
    if target_color==1
        interuptData.targetCloud(t).color = stimulus.colors{1};
        interuptData.otherCloud(t).color = stimulus.colors{2};

    else
        interuptData.targetCloud(t).color = stimulus.colors{2};
        interuptData.otherCloud(t).color = stimulus.colors{1};
    end
    
    if interuptData.nClouds(t)==1  % ignore second cloud
        
        interuptData.targetCloud(t).nDots = randsample(stimulus.minDots:stimulus.max1CloudDots, 1);
        
        interuptData.otherCloud(t).bins = [];
        interuptData.otherCloud(t).nDots = 0;

        [x, y, sizes, dotColors, targetCloud, ~] = gen_dot_info(stimulus, interuptData.targetCloud(t), nan);
        interuptData.otherCloud(t).idx = [];
        
        trialCodeA = 10;
        if length(interuptData.targetCloud(t).bins)==1
            trialCodeB = 1;
        else
            trialCodeB = 3;
        end
        
    else
        interuptData.targetCloud(t).nDots = randsample(stimulus.minDots:stimulus.max2CloudDots, 1);
        
        interuptData.otherCloud(t).bins = get_bins_from_loc_width( ...
            randsample(1:8, 1), ...
            randsample([1, 3], 1), ...
            stimulus.nBins);
        interuptData.otherCloud(t).nDots = randsample(stimulus.minDots:stimulus.max2CloudDots, 1);

        [x, y, sizes, dotColors, targetCloud, otherCloud] = gen_dot_info(stimulus, interuptData.targetCloud(t), interuptData.otherCloud(t));
        interuptData.otherCloud(t).idx = otherCloud.idx; 
        
        trialCodeA = 20;
        if (length(interuptData.targetCloud(t).bins)==1) && (length(interuptData.otherCloud(t).bins)==1)
            trialCodeB = 1;
        elseif (length(interuptData.targetCloud(t).bins)==3) && (length(interuptData.otherCloud(t).bins)==3)
            trialCodeB = 3;
        else
            trialCodeB = 2;
        end
    end

    interuptData.targetCloud(t).idx = targetCloud.idx;
    keep_idx = cat(1, targetCloud.idx, interuptData.otherCloud(t).idx);

    interuptData.disp(t).dotCoords = [x(keep_idx); y(keep_idx)];
    interuptData.disp(t).sizes = sizes(keep_idx);
    interuptData.disp(t).dotColors = dotColors(:, keep_idx);
    
    
    % Interuption info
    interupt_bool = rand() <= stimulus.pInterupt;
    interuptData.interuptCloud(t).bool = interupt_bool;

    interuptData.interuptCloud(t).bins = get_bins_from_loc_width( ...
                randsample(1:8, 1), ... % loc
                1, ... % width==1
                stimulus.nBins);
    interuptData.interuptCloud(t).color = 'grey';
    interuptData.interuptCloud(t).nDots = randsample(stimulus.minDots:stimulus.max1CloudDots, 1); 

    [xInterupt, yInterupt, sizesInterupt, dotColorsInterupt, interuptCloud, ~] = gen_dot_info(stimulus, interuptData.interuptCloud(t), nan);
    interuptData.interuptCloud(t).idx = interuptCloud.idx;
    interuptData.interuptDisp(t).dotCoords = [xInterupt(interuptCloud.idx); yInterupt(interuptCloud.idx)];
    interuptData.interuptDisp(t).sizes = sizesInterupt(interuptCloud.idx);
    interuptData.interuptDisp(t).dotColors = dotColorsInterupt(:, interuptCloud.idx);
        

    
    change = randsample([true, false], 1);
    if change
        trialCodeC = 1;
        
        interuptData.responses(t).CorrectResponse = KbName(stimulus.keyCodes(1));
        probeCloud.bins = interuptData.targetCloud(t).bins;
        probeCloud.color = interuptData.targetCloud(t).color;
        probeCloud.nDots = interuptData.targetCloud(t).nDots;
        probeCloud.idx = targetCloud.idx;

        change_type = randsample([1, 2], 1);
        probeCloud.change_type = change_type;

        % 1 - shift location by 1 bin
        if change_type==1
            shift = randsample([1, -1], 1);

            angles = atan2d(y, x) * pi/180; % retrieve and convert to radians
            angles = angles + shift*(2*pi)/stimulus.nBins; % shift by 1 bin
            dists = sqrt(x.^2 + y.^2);

            x = cos(angles).*dists;
            y = sin(angles).*dists;

        end

        % 2 - change width between 1 and 3
        if change_type==2
            old_bins = probeCloud.bins;
            if length(old_bins)==3
                bins = old_bins(2);
            else
                bins = [old_bins - 1, old_bins, old_bins+1];
            end
            bins = mod(bins, stimulus.nBins);
            bins(bins==0) = stimulus.nBins;
            probeCloud.bins = bins;

            [x, y, sizes, dotColors, newCloud, ~] = gen_dot_info(stimulus, probeCloud, nan);
            probeCloud.idx = newCloud.idx;
        end

    else
        trialCodeC = 0;
        probeCloud = interuptData.targetCloud(t);
        probeCloud.change_type = 0;
        interuptData.responses(t).CorrectResponse = KbName(stimulus.keyCodes(2));
    end
    interuptData.portCode(t) = 200 + 30 + trialCodeA + trialCodeB;

    probeInfo.portCode(t) = 170 + trialCodeC;
    probeInfo.disp(t).dotCoords = [x(probeCloud.idx); y(probeCloud.idx)];
    probeInfo.disp(t).sizes = sizes(probeCloud.idx);
    probeInfo.disp(t).dotColors = dotColors(:, probeCloud.idx); 
    probeInfo.cloud(t) = probeCloud;
end


%% - PTB, Eyetracking Init.

% Psychtoolbox Set-Up
AssertOpenGL;

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');

[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.baseColor);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);

Screen('BlendFunction', ptbWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

[xCenter, yCenter] = RectCenter(winRect);
IFI = Screen('GetFlipInterval', ptbWindow);

[width, height]=Screen('WindowSize', ptbWindow);
stimulus.diodeRect = [width-40, 10,width-10, 40]; %%%% Upper right Corner!!!!!

% PTB Text set-up
Screen('TextFont',ptbWindow,text.instructionFont);
Screen('TextSize',ptbWindow,text.instructionPoints);
Screen('TextStyle',ptbWindow,text.instructionStyle);
ListenChar;

global psych_default_colormode;
psych_default_colormode = 1;

global ptb_drawformattedtext_disableClipping;
ptb_drawformattedtext_disableClipping = 1;

% Introduction Text
RestrictKeysForKbCheck(stimulus.keyCodes);
HideCursor();
clear KbCheck;
introductionText = ['Thank you for participating in this study!\n\n'...
    'In this task, you will see 1 or 2 sets of colored dots - 1 blue and/or 1 green - spread around a central dot.\n'...
    'Afer a delay, one set will reappear.\n\n'...
    'If the set has changed (e.g. in the location or spread of dots), press the '  KbName(stimulus.keyCodes(1)) ' key.\n\n'...
    'During the delay, a set of gray dots will appear briefly. Please ignore the grey dots.'...
    'Press either key to begin.'];
DrawFormattedText(ptbWindow,introductionText,'center','center',text.color,text.instructionWrap);
Screen('Flip',ptbWindow);
[~, ~] = KbWait(-1, 2);


% INITIATE EYE TRACKING
if eyes.tracking
    if EyelinkInit()~= 1;return;end
    
    EyeLinkDefaults=EyelinkInitDefaults(ptbWindow);
    % force remote mode
    if eyes.mode
        Eyelink('command', 'elcl_select_configuration = RTABLER'); % remote mode
        Eyelink('command', 'calibration_type = HV9'); % 9-pt calibration
    else
        Eyelink('command', 'elcl_select_configuration = MTABLER'); % chin rest
        Eyelink('command', 'calibration_type = HV5'); % 5-pt calibration
    end

    Eyelink('command','sample_rate = %d',1000)% set sampling rate
    % stamp header in EDF file
    Eyelink('command', 'add_file_preamble_text','For Spatial Bias Project');
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    screenSize=get(0,'ScreenSize');
    width=screenSize(3);
    height=screenSize(4);
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);
    % make sure that we get gaze data from the Eyelink
    %Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
    % set EDF file contents using the file_sample_data and
    % file-event_filter commands
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERES,STATUS,INPUT');
    % set link data thtough link_sample_data and link_event_filter
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
    % proportion commands to adjust size of calibrated area
    Eyelink('command', 'calibration_area_proportion 0.5 0.5')
    Eyelink('command', 'validation_area_proportion 0.5 0.5')

    % Get host tracker version
    [v,vs]=Eyelink('GetTrackerVersion');

    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    fprintf('Running experiment on version ''%d''.\n', v );

    %Open file to record data to
    eyes.filename=[num2str(participant.num),'_DCD.edf']; %File name MUST BE LESS THAN 8 characters!!!!!
    Eyelink('Openfile', eyes.filename);
end


%% - TEST

RestrictKeysForKbCheck(stimulus.keyCodes);
testIntroductionText = ['Welcome to the main task!.\n\n' ...
    'Remember:\n\n'...
    'If the set has changed, press the '  KbName(stimulus.keyCodes(1)) ' key.\n\n'...
    'If the set has not changed, press the ' KbName(stimulus.keyCodes(2)) ' key.\n\n'...
    'Ignore the set of gray dots that appears during the delay.'...
    'When you are ready to continue, press either key.'];
DrawFormattedText(ptbWindow,testIntroductionText,'center','center',text.color,text.instructionWrap);
Screen('Flip',ptbWindow);
[~, ~] = KbWait(-1, 2);

if eyes.tracking
    RestrictKeysForKbCheck([]);
    EyelinkDoTrackerSetup(EyeLinkDefaults);
    RestrictKeysForKbCheck(stimulus.keyCodes);
end

if porting
    outp(portCode,190);
end
if eyes.tracking
    Eyelink('message','SYNC 190');
end
WaitSecs(.01); % wait ~10msec to separate port codes

breakout = false;
for bN = 1:stimulus.nBlocks    
    if porting
        WaitSecs(.01); % wait 10msec
        outp(portCode,100+bN);
    end
    blockAcc = int16.empty(stimulus.nTrialsPerBlock,0);
    
    for tN = 1:stimulus.nTrialsPerBlock
        noResp = true;
        
        % JITTERED ITI FIXATION
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawingFinished',ptbWindow);
        ITIstart = Screen('Flip',ptbWindow); % show fixation ITI immediately
        if porting
            outp(portCode,150);
        end
        if eyes.tracking
            % set idle mode before start of recording
            Eyelink('Command', 'set_idle_mode');
            % Put block, trial number at the bottom of operater display
            Eyelink('command', 'record_status_message ''BLOCK %d TRIAL %d''', bN, tN);
            % start recording eye position
            % record a few samples before we actually start displaying
            Eyelink('StartRecording'); 
            % Eyelink('message', 'TRIALID %d ', BTRIAL);
            Eyelink('message', 'BLOCK %d ', bN);
            Eyelink('message', 'TRIAL %d ', tN);

            Eyelink('message', 'TrialStart');
        end
        
        % Allow for interruption during the ITI
        RestrictKeysForKbCheck([KbName('e'), KbName('q'), KbName('ESCAPE'), KbName('Return'), KbName('c'), KbName('v')]);
        
        while (GetSecs() - ITIstart) < interuptData.ITIs(tN)
            [keyIsDown, keyTime, keyCode] = KbCheck;
            if keyIsDown
                if strcmp(KbName(keyCode), 'e')
                    EyelinkDoTrackerSetup(EyeLinkDefaults);
                    
                    Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
                    Screen('DrawingFinished',ptbWindow); 
                    Screen('Flip',ptbWindow, GetSecs()+timing.ITIinterruptionWait - (.5*IFI)); % wait 1 sec
                    
                elseif strcmp(KbName(keyCode), 'q')
                    clear KbCheck;
                    RestrictKeysForKbCheck([KbName('y'), KbName('n')]);
                    fprintf('Are you sure you want to end the session? y/n')
                    [keyIsDown, keyCode] = KbWait(-1, 2);
                    if strcmp(KbName(keyCode), 'y')
                        breakout = true;
                        ShowCursor();
                        break;
                    else
                        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
                        Screen('DrawingFinished',ptbWindow); 
                        Screen('Flip',ptbWindow, GetSecs()+timing.ITIinterruptionWait - (.5*IFI)); % wait 1 sec
                    end
                end
            end
            
        end
        clear KbCheck;
        RestrictKeysForKbCheck(stimulus.keyCodes);
        if breakout
            break;
        end
        
        % BASELINE
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawingFinished',ptbWindow);
        baselineStart = Screen('Flip',ptbWindow); % show baseline immediately after jitter is complete
        if porting
            outp(portCode,tN);
        end
        
        if eyes.tracking
            Eyelink('message', 'Baseline');
        end
        
        % MEMORY ARRAY
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawDots', ptbWindow, interuptData.disp(tN).dotCoords, interuptData.disp(tN).sizes, interuptData.disp(tN).dotColors, ...
            [xCenter, yCenter], 2);
        Screen('FillOval',ptbWindow,1,stimulus.diodeRect);
        Screen('DrawingFinished',ptbWindow);
        MemArrayStart = Screen('Flip',ptbWindow,baselineStart+timing.baseline-(.5*IFI));  % show mem array after baseline complete
        if porting
            outp(portCode, interuptData.portCode(tN)) % trial onset, label condition
        end
        if eyes.tracking
            Eyelink('message', ['SYNC ', num2str(interuptData.portCode(tN))]);
        end
        
        
        % DELAY - Pt 1 - 500ms
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawingFinished',ptbWindow);
        Delay1Start = Screen('Flip', ptbWindow, MemArrayStart+timing.memoryArray-(.5*IFI));  % show delay after mem array complete
        if porting
            outp(portCode,160);
        end
        if eyes.tracking
            Eyelink('message','Delay');
        end
        
        % DELAY/Interruption - 150ms
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        if interuptData.interuptCloud(tN).bool
            Screen('DrawDots', ptbWindow, interuptData.interuptDisp(tN).dotCoords, interuptData.interuptDisp(tN).sizes, interuptData.interuptDisp(tN).dotColors, ...
                [xCenter, yCenter], 2);
            PC = 161;
        else
            PC = 162;
        end
        Screen('DrawingFinished',ptbWindow);
        InteruptStart = Screen('Flip', ptbWindow, Delay1Start+timing.postMemoryDelay1-(.5*IFI));  % show interupt/fix after 1st delay period
        if porting
            outp(portCode,PC);
        end
        if eyes.tracking
            Eyelink('message','Delay');
        end
        
        % DELAY - PT 2 - 850ms
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawingFinished',ptbWindow);
        Delay2Start = Screen('Flip', ptbWindow, InteruptStart+timing.postMemoryDelayInterupt-(.5*IFI));  % show 2nd delay after interuption complete
        if porting
            outp(portCode,163);
        end
        if eyes.tracking
            Eyelink('message','Delay');
        end
        
        
        % PROBE & RESPONSE
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawDots', ptbWindow, probeInfo.disp(tN).dotCoords, probeInfo.disp(tN).sizes, probeInfo.disp(tN).dotColors, ...
            [xCenter, yCenter], 2);
        Screen('DrawingFinished',ptbWindow);
        ProbeStart = Screen('Flip', ptbWindow, Delay2Start+timing.postMemoryDelay2-(.5*IFI));  % show probe after 2nd delay period complete
        if porting
            outp(portCode,probeInfo.portCode(tN));
        end
        if eyes.tracking
            Eyelink('message','ProbeDisplay');
        end

        while  noResp
            [keyIsDown, keyTime, keyCode] = KbCheck;
            if keyIsDown
                if porting
                    outp(portCode,180);
                    % draw fixation and wait the minimal amount of time
                    Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
                    Screen('DrawingFinished',ptbWindow); 
                    Screen('Flip',ptbWindow, .005 - (.5*IFI)); % wait 5ms
                end
                if eyes.tracking
                    Eyelink('message','ResponseMade');
                    Eyelink('StopRecording');
                end
                noResp = false;
                interuptData.responses(tN).Response = KbName(keyCode);
                interuptData.responses(tN).ResponseTime = keyTime - ProbeStart;
                blockAcc(tN) = strcmp(interuptData.responses(tN).Response,interuptData.responses(tN).CorrectResponse);
                interuptData.responses(tN).ChoiceAcc = blockAcc(tN);
            end
        end
        fprintf(['b: ' num2str(bN) ' t: ' num2str(tN) ' acc: ' num2str(interuptData.responses(tN).ChoiceAcc) ' RT: ' num2str(interuptData.responses(tN).ResponseTime) '\n']);
        clear KbCheck;
        
    end
    
    if breakout
        break
    end
    
    % BLOCK BREAK
    if bN ~= stimulus.nBlocks
        
        % DISPLAY FEEDBACK
        endBlockText = ['Please take a short break.\n' ...
            'You have completed ' num2str(bN) ' out of ' num2str(stimulus.nBlocks) ' blocks.\n\n'...
            'Your accuracy this block was: ' sprintf('%.2f%%',  mean(blockAcc)*100) '\n\n'...
            'Remember:\n'...
            'If the set has changed, press the '  KbName(stimulus.keyCodes(1)) ' key.\n\n'...
            'If the set has not changed, press the ' KbName(stimulus.keyCodes(2)) ' key.\n\n'...
            'When you are ready to begin the next block, press either key.'];
        DrawFormattedText(ptbWindow, endBlockText, 'center', 'center', text.color, text.instructionWrap);
        Screen('Flip', ptbWindow);

        [~, ~] = KbWait(-1, 2);
        clear KbCheck;
        
        % SAVE DATA SO FAR
        oldDir = cd(paths.behavDir);

        targetCloudTable = struct2table(interuptData.targetCloud);
        targetCloudTable = targetCloudTable(:,ismember(targetCloudTable.Properties.VariableNames, {'bins', 'color', 'nDots'})); % ignore IDX
        targetCloudTable.Properties.VariableNames = {'targetCloud_bins', 'targetCloud_color', 'targetCloud_nDots'};

        interuptCloudTable = struct2table(interuptData.interuptCloud);
        interuptCloudTable = interuptCloudTable(:,ismember(interuptCloudTable.Properties.VariableNames, {'bool', 'bins', 'color', 'nDots'})); % ignore IDX
        interuptCloudTable.Properties.VariableNames = {'interuptCloud_presentBool', 'interuptCloud_bins', 'interuptCloud_color', 'interuptCloud_nDots'};

        otherCloudTable = struct2table(interuptData.otherCloud);
        otherCloudTable = otherCloudTable(:,ismember(otherCloudTable.Properties.VariableNames, {'bins', 'color', 'nDots'})); % ignore IDX
        otherCloudTable.Properties.VariableNames = {'otherCloud_color', 'otherCloud_bins', 'otherCloud_nDots'};

        probeTable = struct2table(probeInfo.cloud);
        probeTable = probeTable(:,ismember(probeTable.Properties.VariableNames, {'bins', 'color', 'nDots', 'change_type'})); % ignore IDX
        probeTable.Properties.VariableNames = {'probe_bins', 'probe_color', 'probe_nDots', 'change_type'};

        portTable = array2table(interuptData.portCode');
        portTable.Properties.VariableNames = {'port_codes'};

        behavTable = [targetCloudTable, interuptCloudTable, otherCloudTable, probeTable, portTable, struct2table(interuptData.responses)];

        writetable(behavTable, [paths.fileBase '.csv']);
        save([paths.fileBase '.mat'], 'participant', 'paths', 'interuptData', 'data', 'probeInfo','stimulus', 'equipment', 'text', 'timing');
        cd(oldDir);
        
        % RECALIBRATE EYE TRACKING
        if eyes.tracking
            RestrictKeysForKbCheck([]);
            EyelinkDoTrackerSetup(EyeLinkDefaults);
            RestrictKeysForKbCheck(stimulus.keyCodes);
        end

        Screen('Flip', ptbWindow);
        WaitSecs(timing.postBreak);
        
    end
end

%% - Saving & Exit

if porting
    WaitSecs(.01); % wait 5=10msec
    outp(portCode,191);
end
if eyes.tracking
    Eyelink('message','SYNC 191');
end
completionText = 'You''ve completed this task! Press either key to end the session.';
DrawFormattedText(ptbWindow, completionText, 'center', 'center', text.color, text.instructionWrap);
Screen('Flip', ptbWindow);
[~, ~] = KbWait(-1, 2);
clear KbCheck;

RestrictKeysForKbCheck([]);
ShowCursor();


% SAVE EYETRACKING
if eyes.tracking
    oldDir = cd(paths.eegDir);
    
    Screen('TextSize',ptbWindow,24);
    Screen('TextFont',ptbWindow,'Arial');
    DrawFormattedText(ptbWindow, 'TRANSFERRING EYE DATA.', 'center', 'center', [255 255 255]);
    Screen('Flip', ptbWindow);
    
    Eyelink('CloseFile');
    
    try
        fprintf('Receiving data file ''%s''\n', eyes.filename );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(eyes.filename, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', eyes.filename, pwd );
        end
    catch rdf
        fprintf('Problem receiving data file ''%s''\n', eyes.filename);
        rdf;
    end
    
    Eyelink('ShutDown');
    cd(oldDir);
end

% % SAVE TASK DATA
oldDir = cd(paths.behavDir);

targetCloudTable = struct2table(interuptData.targetCloud);
targetCloudTable = targetCloudTable(:,ismember(targetCloudTable.Properties.VariableNames, {'bins', 'color', 'nDots'})); % ignore IDX
targetCloudTable.Properties.VariableNames = {'targetCloud_bins', 'targetCloud_color', 'targetCloud_nDots'};

interuptCloudTable = struct2table(interuptData.interuptCloud);
interuptCloudTable = interuptCloudTable(:,ismember(interuptCloudTable.Properties.VariableNames, {'bool', 'bins', 'color', 'nDots'})); % ignore IDX
interuptCloudTable.Properties.VariableNames = {'interuptCloud_presentBool', 'interuptCloud_bins', 'interuptCloud_color', 'interuptCloud_nDots'};

otherCloudTable = struct2table(interuptData.otherCloud);
otherCloudTable = otherCloudTable(:,ismember(otherCloudTable.Properties.VariableNames, {'bins', 'color', 'nDots'})); % ignore IDX
otherCloudTable.Properties.VariableNames = {'otherCloud_color', 'otherCloud_bins', 'otherCloud_nDots'};

probeTable = struct2table(probeInfo.cloud);
probeTable = probeTable(:,ismember(probeTable.Properties.VariableNames, {'bins', 'color', 'nDots', 'change_type'})); % ignore IDX
probeTable.Properties.VariableNames = {'probe_bins', 'probe_color', 'probe_nDots', 'change_type'};

portTable = array2table(interuptData.portCode');
portTable.Properties.VariableNames = {'port_codes'};

behavTable = [targetCloudTable, interuptCloudTable, otherCloudTable, probeTable, portTable, struct2table(interuptData.responses)];

writetable(behavTable, [paths.fileBase '.csv']);
save([paths.fileBase '.mat'], 'participant', 'paths', 'interuptData', 'probeInfo','stimulus', 'equipment', 'text', 'timing');
cd(oldDir);

% Close
sca;
% clear;