clear;
close all;
Screen('CloseAll');
KbName('UnifyKeyNames');

%%  - EEG Session Vars
% screen - 1 for testing room
screenID = 1;

% EEG porting
porting = true;
portCode = hex2dec('D050');

% eyetracking
eyes.tracking = true;
eyes.mode = 0; % 0 = chin rest, 1 = remote

% Practice
do_practice = true;

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
paths.fileBase = ['p-' num2str(participant.num) '_task-changeDetection_Exp1'];

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
timing.postMemoryDelay = 1;
timing.postBreak = 1; 
timing.practiceFeedback = 1;
timing.ITIinterruptionWait = 1;

% Stimulus Details
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

stimulus.nTrialsPractice = 20;
stimulus.nTrials = 1536;
stimulus.nBlocks = 16;
stimulus.nFullLoops = 3;
stimulus.nTrialsPerBlock = stimulus.nTrials/stimulus.nBlocks;

stimulus.sliceThickness = 2*pi / stimulus.nSlices; % radians per slice
stimulus.nDotsPerBin = stimulus.nSlices*stimulus.nRows/stimulus.nBins;
stimulus.binStarts = 1:stimulus.nDotsPerBin:stimulus.nDotsPerBin*(stimulus.nBins-1) + 1;

stimulus.colors = {'blue', 'green'};
stimulus.dotColors.blue =    [0, 0, 255]/255;
stimulus.dotColors.green =   [0, 60, 0]/255;

stimulus.placeholderColor = equipment.baseColor;
stimulus.fixColor = 1;  

stimulus.keyCodes = randsample([KbName('LeftArrow'), KbName('RightArrow')], 2);

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
    %%% 1X = 1 cloud
    %%% 2X = 2 clouds
    %%% X1 = cloud(s) 1 bin width
    %%% X2 = clouds different widths
    %%% X3 = cloud(s) 3 bin width
    
    % Trial section codes
    
    %%% 150 (ITI starts)
    
    %%% trial code (Trial initiated - baseline starts)
    %%% condition code (Stim Display)
   
    %%% 160 (Delay Starts)
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

% PRACTICE TRIALS
% Generated totally at random
practiceData.rng = rng;
practiceData.ITIs = datasample(timing.minITI:timing.stepITI:timing.maxITI, stimulus.nTrials);
for t = 1:stimulus.nTrialsPractice
    practiceData.targetCloud(t).bins = get_bins_from_loc_width( ...
        randsample(1:8, 1), ... % loc
        randsample([1, 3], 1), ... % width
        stimulus.nBins); 
    
    target_color = randsample([1, 2], 1);
    if target_color==1
        practiceData.targetCloud(t).color = stimulus.colors{1};
        practiceData.otherCloud(t).color = stimulus.colors{2};

    else
        practiceData.targetCloud(t).color = stimulus.colors{2};
        practiceData.otherCloud(t).color = stimulus.colors{1};
    end
    
    n_clouds = randsample([1, 2], 1);
    if n_clouds==1  % ignore second cloud
        practiceData.targetCloud(t).nDots = randsample(stimulus.minDots:stimulus.max1CloudDots, 1);
        
        practiceData.otherCloud(t).bins = [];
        practiceData.otherCloud(t).nDots = 0;

        [x, y, sizes, dotColors, targetCloud, ~] = gen_dot_info(stimulus, practiceData.targetCloud(t), nan);
        practiceData.otherCloud(t).idx = [];
    else
        practiceData.targetCloud(t).nDots = randsample(stimulus.minDots:stimulus.max2CloudDots, 1);
        
        practiceData.otherCloud(t).bins = get_bins_from_loc_width( ...
            randsample(1:8, 1), ...
            randsample([1, 3], 1), ...
            stimulus.nBins);
        practiceData.otherCloud(t).nDots = randsample(stimulus.minDots:stimulus.max2CloudDots, 1);

        [x, y, sizes, dotColors, targetCloud, otherCloud] = gen_dot_info(stimulus, practiceData.targetCloud(t), practiceData.otherCloud(t));
        practiceData.otherCloud(t).idx = otherCloud.idx; 
    end

    practiceData.targetCloud(t).idx = targetCloud.idx;
    keep_idx = cat(1, targetCloud.idx, practiceData.otherCloud(t).idx);

    practiceData.disp(t).dotCoords = [x(keep_idx); y(keep_idx)];
    practiceData.disp(t).sizes = sizes(keep_idx);
    practiceData.disp(t).dotColors = dotColors(:, keep_idx);
    
    
    change = randsample([true, false], 1);
    if change
        practiceData.responses(t).CorrectResponse = KbName(stimulus.keyCodes(1));
        probeCloud.bins = practiceData.targetCloud(t).bins;
        probeCloud.color = practiceData.targetCloud(t).color;
        probeCloud.nDots = practiceData.targetCloud(t).nDots;
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
        probeCloud = practiceData.targetCloud(t);
        probeCloud.change_type = 0;
        practiceData.responses(t).CorrectResponse = KbName(stimulus.keyCodes(2));
    end

    practiceProbeInfo.disp(t).dotCoords = [x(probeCloud.idx); y(probeCloud.idx)];
    practiceProbeInfo.disp(t).sizes = sizes(probeCloud.idx);
    practiceProbeInfo.disp(t).dotColors = dotColors(:, probeCloud.idx); 
    practiceProbeInfo.cloud(t) = probeCloud;
end


% TEST TRIALS
data.rng = rng;
data.ITIs = datasample(timing.minITI:timing.stepITI:timing.maxITI, stimulus.nTrials);
image_i = 1;
for l = 1:stimulus.nFullLoops
    for loc_c1 = 1:stimulus.nBins
        for loc_c2 = 1:stimulus.nBins
            for width_c1 = stimulus.binSizes
                for width_c2 = stimulus.binSizes
                    for n_clouds = [1, 2]
                        target_color = randsample([1, 2], 1);
                        change = randsample([true, false], 1);

                        data.targetCloud(image_i).bins = get_bins_from_loc_width(loc_c1, width_c1, stimulus.nBins);
                        

                        if target_color==1
                            data.targetCloud(image_i).color = stimulus.colors{1};
                            data.otherCloud(image_i).color = stimulus.colors{2};

                        else
                            data.targetCloud(image_i).color = stimulus.colors{2};
                            data.otherCloud(image_i).color = stimulus.colors{1};
                        end

                        if n_clouds==1  % ignore second cloud
                            data.targetCloud(image_i).nDots = randsample(stimulus.minDots:stimulus.max1CloudDots, 1);
                            
                            data.otherCloud(image_i).bins = [];
                            data.otherCloud(image_i).nDots = 0;
                            
                            [x, y, sizes, dotColors, targetCloud, ~] = gen_dot_info(stimulus, data.targetCloud(image_i), nan);
                            data.otherCloud(image_i).idx = [];

                            trialCodeA = 10;
                            if width_c1==1
                                trialCodeB = 1;
                            else
                                trialCodeB = 3;
                            end

                        else
                            data.targetCloud(image_i).nDots = randsample(stimulus.minDots:stimulus.max2CloudDots, 1);
                            
                            data.otherCloud(image_i).bins = get_bins_from_loc_width(loc_c2, width_c2, stimulus.nBins);
                            data.otherCloud(image_i).nDots = randsample(stimulus.minDots:stimulus.max2CloudDots, 1);
                            
                            [x, y, sizes, dotColors, targetCloud, otherCloud] = gen_dot_info(stimulus, data.targetCloud(image_i), data.otherCloud(image_i));
                            data.otherCloud(image_i).idx = otherCloud.idx;

                            trialCodeA = 20;
                            if (width_c1==1) && (width_c2==1)
                                trialCodeB = 1;
                            elseif (width_c1==3) && (width_c2==3)
                                trialCodeB = 3;
                            else
                                trialCodeB = 2;
                            end
                        end
                        
                        data.targetCloud(image_i).idx = targetCloud.idx;
                        keep_idx = cat(1, targetCloud.idx, data.otherCloud(image_i).idx);
                        
                        data.disp(image_i).dotCoords = [x(keep_idx); y(keep_idx)];
                        data.disp(image_i).sizes = sizes(keep_idx);
                        data.disp(image_i).dotColors = dotColors(:, keep_idx);
                        
                        if change
                            trialCodeC = 1;
                            data.responses(image_i).CorrectResponse = KbName(stimulus.keyCodes(1));
                            probeCloud.bins = data.targetCloud(image_i).bins;
                            probeCloud.nDots = data.targetCloud(image_i).nDots;
                            probeCloud.color = data.targetCloud(image_i).color;
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
                            probeCloud = data.targetCloud(image_i);
                            probeCloud.change_type = 0;
                            data.responses(image_i).CorrectResponse = KbName(stimulus.keyCodes(2));
                        end
                        data.portCode(image_i) = 200 + trialCodeA + trialCodeB;
                        
                        probeInfo.portCode(image_i) = 170 + trialCodeC;
                        probeInfo.disp(image_i).dotCoords = [x(probeCloud.idx); y(probeCloud.idx)];
                        probeInfo.disp(image_i).sizes = sizes(probeCloud.idx);
                        probeInfo.disp(image_i).dotColors = dotColors(:, probeCloud.idx); 
                        probeInfo.cloud(image_i) = probeCloud;

                        image_i = image_i + 1;
                    end
                end
            end
        end
    end
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
    'If the set has not changed, press the ' KbName(stimulus.keyCodes(2)) ' key.\n\n'...
    'We will begin with a brief practice. Press either key to begin.'];
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

%% - PRACTICE

if do_practice
    
    if eyes.tracking
        RestrictKeysForKbCheck([]);
        EyelinkDoTrackerSetup(EyeLinkDefaults);
        RestrictKeysForKbCheck(stimulus.keyCodes);
    end

    for trial_idx = 1:stimulus.nTrialsPractice
        noResp = true;

        % ITI FIXATION
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        ITIstart = Screen('Flip',ptbWindow); % show fixation ITI immediately

        if eyes.tracking
            % set idle mode before start of recording
            Eyelink('Command', 'set_idle_mode');
            % Put block, trial number at the bottom of operater display
            Eyelink('command', 'record_status_message ''BLOCK %s TRIAL %d''', 'practice', trial_idx);
            % start recording eye position
            % record a few samples before we actually start displaying
            Eyelink('StartRecording'); 
            % Eyelink('message', 'TRIALID %d ', BTRIAL);
            Eyelink('message', 'BLOCK %s ', 'practice');
            Eyelink('message', 'TRIAL %d ', trial_idx);

            Eyelink('Message', 'TrialStart');
        end

        % MEMORY ARRAY
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawDots', ptbWindow, practiceData.disp(trial_idx).dotCoords, practiceData.disp(trial_idx).sizes, practiceData.disp(trial_idx).dotColors, ...
            [xCenter, yCenter], 2);
        MemArrayStart = Screen('Flip',ptbWindow,ITIstart+practiceData.ITIs(trial_idx)-(.5*IFI));  % show mem array after ITI complete
        if eyes.tracking
            Eyelink('message','MemArray');
        end

        % DELAY
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        DelayStart = Screen('Flip', ptbWindow, MemArrayStart+timing.memoryArray-(.5*IFI));  % show delay after mem array complete
        if eyes.tracking
            Eyelink('message','Delay');
        end

        % PROBE & RESPONSE
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawDots', ptbWindow, practiceProbeInfo.disp(trial_idx).dotCoords, practiceProbeInfo.disp(trial_idx).sizes, practiceProbeInfo.disp(trial_idx).dotColors, ...
            [xCenter, yCenter], 2);
        ProbeStart = Screen('Flip', ptbWindow, DelayStart+timing.postMemoryDelay-(.5*IFI));  % show probe after delay period complete
        if eyes.tracking
            Eyelink('message','ProbeDisplay');
        end

        while  noResp
            [keyIsDown, keyTime, keyCode] = KbCheck;
            if keyIsDown
                if porting
                    outp(portCode,14);
                end
                if eyes.tracking
                    Eyelink('message','ResponseMade');
                    Eyelink('StopRecording');
                end
                noResp = false;
                practiceData.responses(trial_idx).Response = KbName(keyCode);
                practiceData.responses(trial_idx).ResponseTime = keyTime - ProbeStart;
            end
        end
        clear KbCheck;

        % FEEDBACK
        if strcmp(practiceData.responses(trial_idx).Response,practiceData.responses(trial_idx).CorrectResponse)
            feedback = 'Correct!';
        else
            feedback = 'Incorrect.';
        end
        DrawFormattedText(ptbWindow,feedback,'center','center',text.color,text.instructionWrap);
        FeedbackStart = Screen('Flip',ptbWindow);
        Screen('Flip', ptbWindow, FeedbackStart+timing.practiceFeedback-(.5*IFI));


    end

    RestrictKeysForKbCheck([KbName('w')]);
    practiceCompletionText = ['Done with practice!\n'...
        'The experimenter will now begin the EEG recording and advance to the main task'];
    DrawFormattedText(ptbWindow,practiceCompletionText,'center','center',text.color,text.instructionWrap);
    Screen('Flip',ptbWindow);
    noCont = true;
    while  noCont
        [keyIsDown, keyTime, keyCode] = KbCheck;
        if keyIsDown
            noCont = false;
        end
    end
    clear KbCheck;
end

%% - TEST

RestrictKeysForKbCheck(stimulus.keyCodes);
testIntroductionText = ['Welcome to the main task!.\n\n' ...
    'Remember:\n\n'...
    'If the set has changed, press the '  KbName(stimulus.keyCodes(1)) ' key.\n\n'...
    'If the set has not changed, press the ' KbName(stimulus.keyCodes(2)) ' key.\n\n'...
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

trial_order = randperm(stimulus.nTrials);
data.trial_order = trial_order;
breakout = false;
for bN = 1:stimulus.nBlocks    
    if porting
        WaitSecs(.01); % wait 10msec
        outp(portCode,100+bN);
    end
    blockAcc = int16.empty(stimulus.nTrialsPerBlock,0);
    
    for tN = 1:stimulus.nTrialsPerBlock
        trial_idx = trial_order(((bN-1)*stimulus.nTrialsPerBlock)+tN);
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
        
        while GetSecs() - ITIstart < data.ITIs(trial_idx)
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
        Screen('DrawDots', ptbWindow, data.disp(trial_idx).dotCoords, data.disp(trial_idx).sizes, data.disp(trial_idx).dotColors, ...
            [xCenter, yCenter], 2);
        Screen('FillOval',ptbWindow,1,stimulus.diodeRect);
        Screen('DrawingFinished',ptbWindow);
        MemArrayStart = Screen('Flip',ptbWindow,baselineStart+timing.baseline-(.5*IFI));  % show mem array after baseline complete
        if porting
            outp(portCode, data.portCode(trial_idx)) % trial onset, label condition
        end
        if eyes.tracking
            Eyelink('message', ['SYNC ', num2str(data.portCode(trial_idx))]);
        end
        
        
        % DELAY
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawingFinished',ptbWindow);
        DelayStart = Screen('Flip', ptbWindow, MemArrayStart+timing.memoryArray-(.5*IFI));  % show delay after mem array complete
        if porting
            outp(portCode,160);
        end
        if eyes.tracking
            Eyelink('message','Delay');
        end
        
        
        % PROBE & RESPONSE
        Screen('DrawDots', ptbWindow, [xCenter, yCenter]', stimulus.fixationSize_pix, stimulus.fixColor, [], 2);
        Screen('DrawDots', ptbWindow, probeInfo.disp(trial_idx).dotCoords, probeInfo.disp(trial_idx).sizes, probeInfo.disp(trial_idx).dotColors, ...
            [xCenter, yCenter], 2);
        Screen('DrawingFinished',ptbWindow);
        ProbeStart = Screen('Flip', ptbWindow, DelayStart+timing.postMemoryDelay-(.5*IFI));  % show probe after delay period complete
        if porting
            outp(portCode,probeInfo.portCode(trial_idx));
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
                data.responses(trial_idx).Response = KbName(keyCode);
                data.responses(trial_idx).ResponseTime = keyTime - ProbeStart;
                blockAcc(tN) = strcmp(data.responses(trial_idx).Response,data.responses(trial_idx).CorrectResponse);
                data.responses(trial_idx).ChoiceAcc = blockAcc(tN);
            end
        end
        fprintf(['b: ' num2str(bN) ' t: ' num2str(tN) ' acc: ' num2str(data.responses(trial_idx).ChoiceAcc) ' RT: ' num2str(data.responses(trial_idx).ResponseTime) '\n']);
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
  
        targetCloudTable = struct2table(data.targetCloud);
        targetCloudTable = targetCloudTable(:,ismember(targetCloudTable.Properties.VariableNames, {'bins', 'color', 'nDots'})); % ignore IDX
        targetCloudTable.Properties.VariableNames = {'targetCloud_bins', 'targetCloud_color', 'targetCloud_nDots'};

        otherCloudTable = struct2table(data.otherCloud);
        otherCloudTable = otherCloudTable(:,ismember(otherCloudTable.Properties.VariableNames, {'bins', 'color', 'nDots'})); % ignore IDX
        otherCloudTable.Properties.VariableNames = {'otherCloud_color', 'otherCloud_bins', 'otherCloud_nDots'};

        probeTable = struct2table(probeInfo.cloud);
        probeTable = probeTable(:,ismember(probeTable.Properties.VariableNames, {'bins', 'color', 'nDots', 'change_type'})); % ignore IDX
        probeTable.Properties.VariableNames = {'probe_bins', 'probe_color', 'probe_nDots', 'change_type'};

        portTable = array2table(data.portCode');
        portTable.Properties.VariableNames = {'port_codes'};

        behavTable = [targetCloudTable, otherCloudTable, probeTable, portTable, struct2table(data.responses)];
        behavTable = behavTable(trial_order, :); % resort by how they actually happened

        writetable(behavTable, [paths.fileBase '.csv']);
        save([paths.fileBase '.mat'], 'participant', 'paths', 'practiceData', 'practiceProbeInfo', 'data', 'probeInfo','stimulus', 'equipment', 'text', 'timing');
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
    WaitSecs(.01); % wait 10msec
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

% SAVE TASK DATA
oldDir = cd(paths.behavDir);
  
targetCloudTable = struct2table(data.targetCloud);
targetCloudTable = targetCloudTable(:,ismember(targetCloudTable.Properties.VariableNames, {'bins', 'color', 'nDots'})); % ignore IDX
targetCloudTable.Properties.VariableNames = {'targetCloud_bins', 'targetCloud_color', 'targetCloud_nDots'};

otherCloudTable = struct2table(data.otherCloud);
otherCloudTable = otherCloudTable(:,ismember(otherCloudTable.Properties.VariableNames, {'bins', 'color', 'nDots'})); % ignore IDX
otherCloudTable.Properties.VariableNames = {'otherCloud_color', 'otherCloud_bins', 'otherCloud_nDots'};

probeTable = struct2table(probeInfo.cloud);
probeTable = probeTable(:,ismember(probeTable.Properties.VariableNames, {'bins', 'color', 'nDots', 'change_type'})); % ignore IDX
probeTable.Properties.VariableNames = {'probe_bins', 'probe_color', 'probe_nDots', 'change_type'};

portTable = array2table(data.portCode');
portTable.Properties.VariableNames = {'port_codes'};

behavTable = [targetCloudTable, otherCloudTable, probeTable, portTable, struct2table(data.responses)];
behavTable = behavTable(trial_order, :); % resort by how they actually happened

writetable(behavTable, [paths.fileBase '.csv']);
save([paths.fileBase '.mat'], 'participant', 'paths', 'practiceData', 'practiceProbeInfo', 'data', 'probeInfo','stimulus', 'equipment', 'text', 'timing');
cd(oldDir);

% Close
sca;
clear;