  function PerceptualGroupingSequentialBalanced

% A task for replicating Irida's sequential task.
% Change detection task with set sizes 2 and 4.
% Each trial contains two stimuli displays.
% Stimuli are sometimes divided over the two displays.
% Stimuli on the second display appear either
% in the same locations used in the first display
% or in different locations.
% Gray circles are placeholder to control for visual
% stimulation across set sizes and location conditions.
% Gray color matched to be the average luminance of 
% color set.

close all;
clc;
dbstop if error % tell us what the error is if there is one
AssertOpenGL; % make sb ure openGL rendering is working (aka psychtoolbox is on the path)
KbName('UnifyKeyNames');
%-------------------------------------------------------------------------
% Build a GUI to get subject number
%-------------------------------------------------------------------------
prompt = {'USI','Subject Number','Age','Sex (Female, Male, Other'}; % what info do we want from the subject?
defAns = {'','','',''};                                             % fill in some default answers
box = inputdlg(prompt,'Enter Subject Information', 1, defAns);      % build the GUI

p.clockOutput = clock;                          % record time and date!! 

if length(box) == length(defAns)                % check for enough input or bail
    p.subUSI = box{1};
    p.subNum = box{2};
    p.subAge = box{3};
    p.subSex = num2str(box{4});
else
    return;
end

%-------------------------------------------------------------------------
% Important options for all experiments
%-------------------------------------------------------------------------
p.is_PC = ispc; % ispc function detects if it's a pc computer or not
p.portCodes = 1; % 1 = use p.portCodes (we're in the booth)
p.windowed = 0; % 1 = small win for easy debugging!
p.eyeTrack = 1; % 1 = eyetracking

if p.eyeTrack
    p.eyeMode = 0; % 1 = remote, 0 = chin rest
    % for remote use XX mm lens; monitor -> eye dist = 100, camera ->
    % eye dist = 60 (below 55, above 62 ish it freaks out)
    
    % CHIN REST MONOCULAR (P. 59 IN MANUAL) use 35 mm lens; monitor -> eye dist = less imiportant (the sub-portion of the monitor
    % that you want to track should subtend no more than 32 X 25 degrees of visual angle, and eye-to-monitor distance,
    % should be at least 1.75 times the effective display width (e.g. width you want to track)
    % camera -> eye dist = 40c - 70 cm. The ideal distance is 50 to 55 cm.
end
%-------------------------------------------------------------------------
% Build an output directory & check to make sure it doesn't already exist
%-------------------------------------------------------------------------
p.root = pwd;
% if the subject data directory doesn't exist, make one!!
if ~exist([p.root,filesep,'Subject Data',filesep], 'dir');
    mkdir([p.root,filesep,'Subject Data',filesep]);
end
%-------------------------------------------------------------------------
% Build psychtoolbox window & hide the task bar
%-------------------------------------------------------------------------
win = openWindow(p);
% Manually hide the task bar so it doesn't pop up
if p.is_PC
    ShowHideWinTaskbarMex(0);
end
%-------------------------------------------------------------------------
% Run Experiment
%-------------------------------------------------------------------------
% Build an output file and check to make sure that it doesn't exist yet
fileName = [p.root,filesep,'Subject Data',filesep,num2str(p.subNum), '_PerceptualGroupingSequentialBalanced.mat'];
if p.subNum ~= 0
    if exist(fileName,'file')
        Screen('CloseAll');
        msgbox('File already exists!', 'modal')
        return;
    end
end
%----------------------------------------------------
% Port Settings - % booth 1 p.portCodes: event DCC8 / response
%----------------------------------------------------
p.portCodeNumbers = [11:15,19:23,100:106];
p.portCodeLabels = {'SS2 Same First','SS2 Same Second','SS2 Different First',...
    'SS2 Different Second','SS2 Simultaneous',...
    'SS4 Same First','SS4 Same Second',...
    'SS4 Different First','SS4 Different Second','SS4 Simultaneous',...
    'ITI','Stimuli 1','Retention 1','Stimuli 2','Retention 2','Test','Response'};

%%%% Guide to memory array  portCodes
% 11 - SS2 Same First
% 12 - SS2 Same Second
% 13 - SS2 Different First
% 14 - SS2 Different Second
% 15 - SS2 Simultaneous
% 19 - SS4 Same First
% 20 - SS4 Same Second
% 21 - SS4 Different First
% 22 - SS4 Different Second
% 23 - SS4 Simultaneous

if p.portCodes == 1
    % run the script to configure the parallel port
    config_io;
    % write a value to the default LPT1 printer output port (at 0x378)
    % This is the port address!
end
event_port = hex2dec('D050');
%----------------------------------------------------
% Get screen params, build the display
%----------------------------------------------------
% commandwindow; % select the command win to avoid typing in open scripts
ListenChar(2); % don't print things in the command window

% Set the random state to the random seed at the beginning of the experiment!!
rng default; rng shuffle; p.rngSettings = rng; % alt to rand with 'state' input
Screen('TextFont', win.onScreen, 'Arial');
prefs = getPreferences(); % function that grabs all of our preferences

% Set up fixation point rect (b/c uses both prefs and win)
win.fixRect = [(win.centerX - prefs.fixationSize),(win.centerY - prefs.fixationSize), ...
    (win.centerX  + prefs.fixationSize), (win.centerY + prefs.fixationSize)];

while KbCheck; end;
KbName('UnifyKeyNames');   % This command switches keyboard mappings to the OSX naming scheme, regardless of computer.
% unify key names so we don't need to mess when switching from mac
% to pc ...
escape = KbName('ESCAPE');  % Mac == 'ESCAPE' % PC == 'esc'
prefs.changeKey = KbName('/?'); % on mac, 56 % 191 == / pc
prefs.nochangeKey = KbName('z'); % on mac, 29  % 90 == z
space = KbName('space');
%--------------------------------------------------------
% Preallocate some variable structures! :)
%--------------------------------------------------------
% Stimulus params
stim.condition = NaN(prefs.numTrials,prefs.numBlocks); % Same or different locations
stim.change = NaN(prefs.numTrials,prefs.numBlocks);
stim.setSize = NaN(prefs.numTrials,prefs.numBlocks); % 2, 4, or 6
stim.probeItem = NaN(prefs.numTrials,prefs.numBlocks);
stim.probeDisplay = NaN(prefs.numTrials,prefs.numBlocks);
stim.colors1 = cell(prefs.numTrials,prefs.numBlocks);
stim.colors2 = cell(prefs.numTrials,prefs.numBlocks);
stim.probeColor = NaN(prefs.numTrials,prefs.numBlocks); % color presented during the actual probe test
% Response vectors
stim.response = NaN(prefs.numTrials,prefs.numBlocks);
stim.accuracy = NaN(prefs.numTrials,prefs.numBlocks);
stim.rt = NaN(prefs.numTrials,prefs.numBlocks);
stim.x1 = cell(prefs.numTrials,prefs.numBlocks);
stim.x2 = cell(prefs.numTrials,prefs.numBlocks);
stim.y1 = cell(prefs.numTrials,prefs.numBlocks);
stim.y2 = cell(prefs.numTrials,prefs.numBlocks);
stim.center = cell(prefs.numTrials,prefs.numBlocks);
%--------------------------------------------------------
% Preallocate space for timing parameters
%--------------------------------------------------------
time.iti.vblstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.iti.onset =  NaN(prefs.numTrials,prefs.numBlocks);
time.iti.flipstamp = NaN(prefs.numTrials,prefs.numBlocks);
time.iti.missed = NaN(prefs.numTrials,prefs.numBlocks);

time.stim1.vblstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.stim1.onset =  NaN(prefs.numTrials,prefs.numBlocks);
time.stim1.flipstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.stim1.missed =  NaN(prefs.numTrials,prefs.numBlocks);

time.retention1.vblstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.retention1.onset = NaN(prefs.numTrials,prefs.numBlocks);
time.retention1.flipstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.retention1.missed = NaN(prefs.numTrials,prefs.numBlocks);

time.stim2.vblstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.stim2.onset =  NaN(prefs.numTrials,prefs.numBlocks);
time.stim2.flipstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.stim2.missed =  NaN(prefs.numTrials,prefs.numBlocks);

time.retention2.vblstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.retention2.onset = NaN(prefs.numTrials,prefs.numBlocks);
time.retention2.flipstamp =  NaN(prefs.numTrials,prefs.numBlocks);
time.retention2.missed = NaN(prefs.numTrials,prefs.numBlocks);

%---------------------------------------------------
% Put up instructions
instruct(win)
%----------------------------------------------------
% Initiate the eyetracker
%----------------------------------------------------
% Initialization of the connection with Eyelink
% Exit program if this fails
if p.eyeTrack
    if EyelinkInit()~= 1;return;end
end
%----------------------------------------------------
% Now that the display is built, send experiment details to the eye tracker
%----------------------------------------------------
if p.eyeTrack
    EyeLinkDefaults=EyelinkInitDefaults(win.onScreen);
    % force remote mode    
    if p.eyeMode
        Eyelink('command', 'elcl_select_configuration = RTABLER'); % remote mode
        Eyelink('command', 'calibration_type = HV5'); % 5-pt calibration  
    else
        Eyelink('command', 'elcl_select_configuration = MTABLER'); % chin rest
        Eyelink('command', 'calibration_type = HV9'); % 9-pt calibration        
    end    
    
    Eyelink('command','sample_rate = %d',1000)% set sampling rate
    % stamp header in EDF file
    Eyelink('command', 'add_file_preamble_text','For Perceptual Grouping Project');
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    %%[width, height]=Screen('WindowSize', s);
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
    edfFile=[num2str(p.subNum),'_PGSB.edf']; %File name MUST BE LESS THAN 8 characters!!!!!
    Eyelink('Openfile', edfFile);
end
tic

%---------------------------------------------------
% Begin Block Loop
%---------------------------------------------------
for b = 1:prefs.numBlocks

    %-----------------------------------------
    % Calibrate the eyetracker
    %-----------------------------------------
%     if mod(b,2) == 1 % after every two blocks
        if p.eyeTrack
            EyelinkDoTrackerSetup(EyeLinkDefaults);
        end
%     end
    
    instruct_block(win)
        
    % Pick out the order of trials for this block, based on full Factorial Design
    prefs.order(:,b) = Shuffle(1:prefs.numTrials);
    stim.condition(:,b) = prefs.conditions(prefs.fullFactorialDesign(prefs.order(:,b),1));
    stim.setSize(:,b) = prefs.setSizes(prefs.fullFactorialDesign(prefs.order(:,b),2));
    stim.probeDisplay(:,b) = prefs.probeDisplay(prefs.fullFactorialDesign(prefs.order(:,b),3));
    stim.change(:,b) = prefs.change(prefs.fullFactorialDesign(prefs.order(:,b),4));
    
    %-------------------------------------------------------
    % Begin Trial Loop
    %-------------------------------------------------------
    for t = 1:prefs.numTrials
        %--------------------------------------------------------
        % Figure out the conditions for this trial!
        %--------------------------------------------------------
        condition = stim.condition(t,b);
        setSize = stim.setSize(t,b);
        if condition == 3 % simultaneous condition only probes first display
            stim.probeDisplay(t,b) = 1;
        end
        probeDisplay = stim.probeDisplay(t,b);
        
        change = stim.change(t,b);
        if condition == 3
            nItems = setSize;
            changeIndex = randperm(setSize);
        else
            nItems = setSize/2;
            changeIndex = randperm(setSize/2);
        end
        stim.probeItem(t,b) = changeIndex(1);
        
        ITI = .01*randi([60 100],1); % random interval of 600 - 1000 ms
        
        change
        probeDisplay
        condition
        setSize
        
        % calculate locations
        [xPos1,yPos1,xPos2,yPos2,cxPos1,cyPos1,cxPos2,cyPos2,firstDispCenters] = getStimLocs(prefs,win,setSize,condition);
        
        stim.x1{t,b} = xPos1; stim.y1{t,b} = yPos1;
        stim.x2{t,b} = xPos2; stim.y2{t,b} = yPos2;
        stim.cx1{t,b} = cxPos1; stim.cy1{t,b} = cyPos1;
        stim.cx2{t,b} = cxPos2; stim.cy2{t,b} = cyPos2;
        stim.center{t,b} = firstDispCenters;
        
        colorIndex1 = randperm(size(win.colors_9,1));
        stim.colors1{t,b} = colorIndex1(1:nItems);
        if condition ~= 3
            colorIndex2 = colorIndex1(nItems+1:end);
            stim.colors2{t,b} = colorIndex2(1:nItems);
        end
        
        if probeDisplay == 1
            sColor = colorIndex1(changeIndex(1));  % sColor is the square-of-interest's color if NOT a change condition!
            dColors = Shuffle(colorIndex1(~ismember(colorIndex1,sColor))); % different colors from chosen square
            
            if condition == 1
                dColors = dColors(~ismember(dColors,colorIndex2(changeIndex(1))));
            end
            
            changeLocX = xPos1(changeIndex(1));
            changeLocY = yPos1(changeIndex(1));
        else
            if condition ~= 3
                sColor = colorIndex2(changeIndex(1));  % sColor is the square-of-interest's color if NOT a change condition!
                dColors = Shuffle(colorIndex2(~ismember(colorIndex2,sColor))); % different colors from chosen square
                changeLocX = xPos2(changeIndex(1));
                changeLocY = yPos2(changeIndex(1));
            end
        end
        changeColor = dColors(1); % now we use the index to pick the change color!
        
        if change == 1
            testColor = changeColor;
        else
            testColor = sColor;
        end
        stim.probeColor(t,b) = testColor;
        
%         stim.changeInd(t,b) = changeIndex;       
        
        %--------------------------------------------------------
        % Create and flip up the ITI display
        %--------------------------------------------------------
        Screen('FillRect',win.onScreen,win.foreColor,win.foreRect);      % Draw the foreground win
        Screen('FillOval',win.onScreen,win.black,win.fixRect);           % Draw the fixation point
        Screen('DrawingFinished',win.onScreen);                          % Tell ptb we're done drawing for the moment (makes subsequent flip command execute faster)

        % Flip ITI
        [time.iti.vblstamp(t,b),...
            time.iti.onset(t,b),...
            time.iti.flipstamp(t,b),...
            time.iti.missed(t,b)] ...
            = Screen('Flip',win.onScreen);
        
        %-----------------------------------------------------------------------
        % Send eyetracking message and EEG trigger for ITI display
        %-----------------------------------------------------------------------
        if p.eyeTrack
            % set idle mode before start of recording
            Eyelink('Command', 'set_idle_mode');
            % Put block, trial number at the bottom of operater display
            Eyelink('command', 'record_status_message ''BLOCK %d TRIAL %d''', b, t);
            % start recording eye position
            % record a few samples before we actually start displaying
            Eyelink('StartRecording'); % KIRSTEN SHOULD MOVE THIS!!!!
            %       Eyelink('message', 'TRIALID %d ', BTRIAL);
            Eyelink('message', 'BLOCK %d ', b);
            Eyelink('message', 'TRIAL %d ', t);

            Eyelink('Message', 'TrialStart');
        end
        
        % EEG trigger for ITI display
        if p.portCodes == 1
            outp(event_port,100);
        end

        % Wait the ITI
        WaitSecs('UntilTime',time.iti.onset(t,b)+ITI-(0.5*prefs.refreshCycle));

        %-----------------------------------------------------------------------
        % Create and flip up the first stimuli display
        %-----------------------------------------------------------------------
        Screen('FillRect',win.onScreen,win.foreColor,win.foreRect);
        Screen('FillOval',win.onScreen,win.black,win.fixRect);
        Screen('FillOval',win.onScreen,win.white,win.probeRect);     % Draw the stimtrak probe win
        Screen('FillRect',win.onScreen,win.colors_9(colorIndex1(1:nItems),:)',[(xPos1-prefs.stimSize/2);...
            (yPos1-prefs.stimSize/2);(xPos1+prefs.stimSize/2);(yPos1+prefs.stimSize/2)]);
        Screen('FillOval',win.onScreen,win.cGrey,[(cxPos1-prefs.circleSize);...
            (cyPos1-prefs.circleSize);(cxPos1+prefs.circleSize);(cyPos1+prefs.circleSize)]);
        Screen('DrawingFinished',win.onScreen);
        
        % Flip stimulus display
        [time.stim1.vblstamp(t,b),...
            time.stim1.onset(t,b),...
            time.stim1.flipstamp(t,b),...
            time.stim1.missed(t,b)] ...
            = Screen('Flip',win.onScreen);

        %-----------------------------------------------------------------------
        % Send EEG trigger and eyetracking message for the first stimulus display
        %-----------------------------------------------------------------------
        if p.portCodes == 1
            outp(event_port,(condition*2)+(setSize*4)+probeDisplay);
        end

        if p.eyeTrack
            Eyelink('message','StimuliArray');
        end

        % Wait the stimulus duration
        WaitSecs('UntilTime',time.stim1.onset(t,b)+prefs.stimulusDuration-(0.5*prefs.refreshCycle));
        
        %-----------------------------------------------------------------------
        % Create and flip up the first retention interval
        %-----------------------------------------------------------------------
        % Draw blank screen
        Screen('FillRect',win.onScreen,win.foreColor,win.screenRect);    % Draw the foreground win
        Screen('FillRect',win.onScreen,win.foreColor,win.foreRect);      % Draw the foreground win
        Screen('FillOval',win.onScreen,win.black,win.fixRect);           % Draw the fixation point
        Screen('DrawingFinished',win.onScreen);

        % Flip blank retention screen
        [time.retention1.vblstamp(t,b),...
            time.retention1.onset(t,b),...
            time.retention1.flipstamp(t,b),...
            time.retention1.missed(t,b)] ...
            = Screen('Flip',win.onScreen);

        %-----------------------------------------------------------------------
        % Send EEG trigger and eyetracking message for the first retention interval
        %-----------------------------------------------------------------------
        if p.portCodes == 1
            outp(event_port,102);
        end

        if p.eyeTrack
            Eyelink('message','Retention');
        end

        % Wait the retention interval
        WaitSecs('UntilTime',time.retention1.onset(t,b)+prefs.interStimInterval-(0.5*prefs.refreshCycle));
        
        %-----------------------------------------------------------------------
        % Create and flip up the second stimuli display
        %-----------------------------------------------------------------------
        Screen('FillRect',win.onScreen,win.foreColor,win.foreRect);
        Screen('FillOval',win.onScreen,win.black,win.fixRect);
        Screen('FillOval',win.onScreen,win.white,win.probeRect);     % Draw the stimtrak probe win
        if condition ~= 3 
            Screen('FillRect',win.onScreen,win.colors_9(colorIndex2(1:nItems),:)',[(xPos2-prefs.stimSize/2);...
                (yPos2-prefs.stimSize/2);(xPos2+prefs.stimSize/2);(yPos2+prefs.stimSize/2)]);
            Screen('FillOval',win.onScreen,win.cGrey,[(cxPos2-prefs.circleSize);...
                (cyPos2-prefs.circleSize);(cxPos2+prefs.circleSize);(cyPos2+prefs.circleSize)]);
        end
        Screen('DrawingFinished',win.onScreen);

        % Flip stimulus display
        [time.stim2.vblstamp(t,b),...
            time.stim2.onset(t,b),...
            time.stim2.flipstamp(t,b),...
            time.stim2.missed(t,b)] ...
            = Screen('Flip',win.onScreen);
        
        %-----------------------------------------------------------------------
        % Send EEG trigger and eyetracking message for the second stimulus display
        %-----------------------------------------------------------------------
        if p.portCodes == 1
            outp(event_port,103);
        end

        if p.eyeTrack
            Eyelink('message','StimuliArray2');
        end

        % Wait the stimulus duration
        WaitSecs('UntilTime',time.stim2.onset(t,b)+prefs.stimulusDuration-(0.5*prefs.refreshCycle));
        
        %-----------------------------------------------------------------------
        % Create and flip up the second retention interval
        %-----------------------------------------------------------------------
        % Draw blank screen
        Screen('FillRect',win.onScreen,win.foreColor,win.screenRect);    % Draw the foreground win
        Screen('FillRect',win.onScreen,win.foreColor,win.foreRect);      % Draw the foreground win
        Screen('FillOval',win.onScreen,win.black,win.fixRect);           % Draw the fixation point
        Screen('DrawingFinished',win.onScreen);

        % Flip blank retention screen
        [time.retention2.vblstamp(t,b),...
            time.retention2.onset(t,b),...
            time.retention2.flipstamp(t,b),...
            time.retention2.missed(t,b)] ...
            = Screen('Flip',win.onScreen);

        %-----------------------------------------------------------------------
        % Send EEG trigger and eyetracking message for the second retention interval
        %-----------------------------------------------------------------------
        if p.portCodes == 1
            outp(event_port,104);
        end

        if p.eyeTrack
            Eyelink('message','Retention2');
        end

        % Wait the retention interval
        WaitSecs('UntilTime',time.retention2.onset(t,b)+prefs.retentionInterval-(0.5*prefs.refreshCycle));
        
        %-----------------------------------------------------------------------
        % Send EEG trigger and eyetracking message for the test display
        %-----------------------------------------------------------------------
        if p.portCodes == 1
            outp(event_port,105);
        end

        if p.eyeTrack
            Eyelink('message', 'TestArray');
        end

        % Wait for a response
        rtStart = GetSecs;        

        Screen('FillRect',win.onScreen,win.foreColor,win.foreRect);
        Screen('FillOval',win.onScreen,win.black,win.fixRect);
        Screen('FillRect',win.onScreen,win.colors_9(testColor,:)',[(changeLocX-prefs.stimSize/2),(changeLocY-prefs.stimSize/2),(changeLocX+prefs.stimSize/2),(changeLocY+prefs.stimSize/2)]);
        Screen('DrawingFinished',win.onScreen);
        Screen('Flip',win.onScreen);
        
        while 1
            [keyIsDown,secs,keyCode]=KbCheck;
            if keyIsDown
                if keyCode(escape)                              % if escape is pressed, bail out
                    ListenChar(0);
                    % save data file if we abort the session
                    p.date_time_end = clock; % record time and date of the end of the sesssion
                    save(fileName,'p','stim','prefs','win','time');
                    Screen('CloseAll');
                    return;
                end
                kp = find(keyCode);
                kp = kp(1); % in case they press 2 buttons at the EXACT same time!!! Not that this happened at the most aggravating possible point in some previous experiment sessions, but yep. it did.
                if kp== prefs.changeKey || kp== prefs.nochangeKey
                    stim.response(t,b)=kp;
                    rtEnd = GetSecs;
                    break
                end
            end
        end
        
        % Blank the screen after the response 
        Screen('FillRect',win.onScreen,win.foreColor);            % Draw the foreground win
        Screen('FillOval',win.onScreen,win.black,win.fixRect);           % Draw the fixation point
        Screen('DrawingFinished',win.onScreen);
        Screen('Flip',win.onScreen);
        
        stim.rt(t,b) = rtEnd-rtStart;

        if change == 1
            if stim.response(t,b) == prefs.changeKey
                stim.accuracy(t,b)=1;
            else
                stim.accuracy(t,b)=0;
            end
        else
            if stim.response(t,b) == prefs.nochangeKey
                stim.accuracy(t,b)=1;
            else
                stim.accuracy(t,b)=0;
            end
        end
        stim.accuracy(t,b)
        
        % EEG trigger for response press
        if p.portCodes == 1
            outp(event_port,106);
        end        
        
        WaitSecs(.1); % wait a little to prevent accidental starting of next trial when spacebar press is long.

        Screen('FillRect',win.onScreen,win.foreColor,win.foreRect);      % Draw the foreground win
        Screen('FillOval',win.onScreen,win.black,win.fixRect);           % Draw the fixation point
        Screen('DrawingFinished',win.onScreen);
        Screen('Flip',win.onScreen);

        %-----------------------------
        % End recording of eyetracker
        %-----------------------------
        if p.eyeTrack
            Eyelink('message', 'Response');
            Eyelink('StopRecording');
        end


        % save data file at the end of each trial
        Screen('Close');

        save(fileName,'p','stim','prefs','win','time');
        
    end % end of trial loop
    
    % Save data file at the end of each block
    save(fileName,'p','stim','prefs','win','time');

    % Tell subjects that they've finished the current block / the experiment
    if b<prefs.numBlocks
        tic
        while toc < prefs.breakLength*60;
            tocInd = round(toc);
            Screen('FillRect',win.onScreen,win.foreColor,win.foreRect);      % Draw the foreground win
            Screen('FillOval',win.onScreen,win.black,win.fixRect);           % Draw the fixation point
            Screen('FillOval',win.onScreen,win.foreColor,win.probeRect);     % Draw the stimtrak probe win
            total_acc = sum(stim.accuracy(:,b));
            accuracy = round((total_acc/prefs.numTrials)*100);
            DrawFormattedText(win.onScreen, sprintf('%d Percent Correct', accuracy), 'center', win.centerY+90, [255 255 255]);
            DrawFormattedText(win.onScreen, 'Take a break.','center',win.centerY-90,win.white);
            DrawFormattedText(win.onScreen, ['Time Remaining: ',char(num2str((prefs.breakLength*60)-tocInd))],'center',win.centerY-60,[255 0 0]);
            DrawFormattedText(win.onScreen, ['Block ',num2str(b),' of ',num2str(prefs.numBlocks),' completed.'],'center',win.centerY+60,win.white);

            Screen('Flip', win.onScreen);
        end
    end

    if b == prefs.numBlocks;

        Screen('TextSize',win.onScreen,24);
        Screen('TextFont',win.onScreen,'Arial');
        DrawFormattedText(win.onScreen, 'Finished! Please wait for the experimenter.', 'center', 'center', [255 255 255]);
        Screen('Flip', win.onScreen);

        while KbCheck; end;
        KbName('UnifyKeyNames');
        space = KbName('space');
        while 1
            [keyIsDown,s,keyCode]=KbCheck;
            if keyIsDown
                kp = find(keyCode);
                if kp == space
                    break;
                end
            end
        end
    end
end    % end of the block loop

%-----------------------------
% Close eyetracking file
%-----------------------------
if p.eyeTrack
    Screen('TextSize',win.onScreen,24);
    Screen('TextFont',win.onScreen,'Arial');
    DrawFormattedText(win.onScreen, 'TRANSFERRING EYE DATA.', 'center', 'center', [255 255 255]);
    Screen('Flip', win.onScreen);
    
    Eyelink('CloseFile');
    
    try
    fprintf('Receiving data file ''%s''\n', edfFile );
    status=Eyelink('ReceiveFile');
    if status > 0
        fprintf('ReceiveFile status %d\n', status);
    end
    if 2==exist(edfFile, 'file')
        fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
    end
    catch rdf
        fprintf('Problem receiving data file ''%s''\n', edfFile );
        rdf;
    end
    
%     cleanup;
    
    Eyelink('ShutDown');
end

%-------------------------------------------------------------------------
% Close psychtoolbox window & Postpare the environment
%-------------------------------------------------------------------------
sca;
ListenChar(0);
if p.is_PC
ShowHideWinTaskbarMex(1);
end
close all;
clear all;
end
%-------------------------------------------------------------------------
function instruct(win)

InstructImage = imread([pwd,'/Instructions_Sequential_Balanced'],'png','BackgroundColor',[win.gray/255,win.gray/255,win.gray/255]);
InstructImage = imresize(InstructImage,0.5);

textOffset = 300;
textSize = 20;
% textFont = 'Arial';

sizeInstruct = size(InstructImage);
rectInstruct = [0 0 sizeInstruct(2) sizeInstruct(1)];
rectTestCoor = [win.centerX,win.centerY-150];


InstructText = ['1. Keep your eyes on the center fixation circle. \n'...
    '2. Remember the color(s) of the square(s). \n'...
    '3. Keep the colors in mind. \n'...
    '4. Remember the color(s) of the squares in the second display if there is one, \n'...
    'which might occupy the same locations as the previous squares. \n'...
    '5. Keep all colors in mind. \n'...
    '5. One of the squares will reappear. Is it the same color as the \n'...
    'square(s) that appeared in that location? \n'...
    '  \n'...
    'Keep your eyes on the center fixation circle for the entire trial. \n'...
    'Blink on the response screen. \n'...
    '  \n'...
    'Press the spacebar to continue.'];

% Show image with explanatory text
Screen('FillRect', win.onScreen, win.gray);
Screen('TextSize', win.onScreen, win.fontsize);

Screen('PutImage',win.onScreen,InstructImage,CenterRectOnPoint(rectInstruct,rectTestCoor(1),rectTestCoor(2)));

Screen('TextSize', win.onScreen, textSize);

DrawFormattedText(win.onScreen, InstructText, win.centerX-textOffset,win.centerY+(sizeInstruct(1)*.35),win.white,[],[],[],1.25);

Screen('Flip', win.onScreen);

% Wait for a spacebar press to continue with next block
while KbCheck; end;
KbName('UnifyKeyNames');   % This command switches keyboard mappings to the OSX naming scheme, regardless of computer.
space = KbName('space');
while 1
    [keyIsDown,s,keyCode]=KbCheck;
    if keyIsDown
        kp = find(keyCode);
        if kp == space
            break;
        end
    end
end

end
%-------------------------------------------------------------------------
function instruct_block(win)

textSize = 25;

% Show image with explanatory text
Screen('FillRect', win.onScreen, win.gray);
Screen('TextSize', win.onScreen, win.fontsize);
Screen('TextSize', win.onScreen, textSize); % 24 = number pixels
DrawFormattedText(win.onScreen, 'Press the spacebar to begin the block.','center','center',win.white);
Screen('Flip', win.onScreen);

% Wait for a spacebar press to continue with the block
while KbCheck; end;
KbName('UnifyKeyNames');
space = KbName('space');
while 1
    [keyIsDown,s,keyCode]=KbCheck;
    if keyIsDown
        kp = find(keyCode);
        if kp == space
            break;
        end
    end
end
end
%-------------------------------------------------------------------------
function [xPos1,yPos1,xPos2,yPos2,cxPos1,cyPos1,cxPos2,cyPos2,firstDispCenters] = getStimLocs(prefs,win,setSize,condition)

clear xPos1 yPos1 xPos2 yPos2 cxPos1 cyPos1 cxPos2 cyPos2

xCenter = win.centerX;
yCenter = win.centerY;

m = 2; % number of columns
n = m; % number of rows

xspacing = 200; 
yspacing = 200;

diag_length = round(sqrt(xspacing^2 + yspacing^2)); %diagonal of grid points
% dis_center = 150;
dis_center = diag_length/2; %distance between fixation point and corner of grid

if mod(m,2)==0
    cDis = (((m/2)-1) * diag_length) + diag_length/2 + dis_center;
else
    cDis = ((m-1)/2 * diag_length) + dis_center;
end

grid_width = (m-1)/2;
grid_height = (n-1)/2;

xCoords = linspace(-grid_width,grid_width,m);
yCoords = linspace(-grid_height,grid_height,n);

cnt = 1;
for yPos = 1:n
    for xPos = 1:m
       x(cnt,1) = xCoords(xPos);
       y(cnt,1) = yCoords(yPos);
       cnt = cnt +1;
    end
end

% scale size of grid
x = xspacing*x;
y = yspacing*y;

qAng = 45:90:360;
for q = 1:length(qAng)
    xGrid(:,q) = round(x + (cDis*cosd(qAng(q))));
    yGrid(:,q) = round(y + (cDis*sind(qAng(q))));
end

xGrid = xGrid + xCenter;
yGrid = yGrid + yCenter;

xmax = max(max(xGrid));
ymax = max(max(yGrid));

xlim = win.screenX - win.screenX * .2;
ylim = win.screenY - win.screenY * .2;

if xmax > xlim || ymax > ylim
    fprintf('Error: grid locations exceed screen limits \n')
    return;
end

nQuads = 1:4; % no. of quadrants
p = m*n*length(nQuads);
nCenters = 1:p; % no. of possible locations in each quadrant, m*n
nStim = 4; %constant total no. of items on each display

% firstDispQuads = datasample(nQuads,setSize/2,'Replace',true); % randomly pick quadrant(s) for first display w/replacement
firstDispCenters = datasample(nCenters,nStim,'Replace',false); % randomly pick stimuli locs
if iscolumn(firstDispCenters)
    firstDispCenters = firstDispCenters'; % transpose if it's a col bc it doesn't work with nCenter = 1
end

xGrid = xGrid(:);
yGrid = yGrid(:);

if condition ~= 3
    for l = 1:setSize/2
        xPos1(l) = xGrid(firstDispCenters(l));
        yPos1(l) = yGrid(firstDispCenters(l));
    end

    for k = 1:(nStim - setSize/2)
        cxPos1(k) = xGrid(firstDispCenters(k+setSize/2));   %x and y pos of placeholder circles
        cyPos1(k) = yGrid(firstDispCenters(k+setSize/2));
    end
end
    
if condition == 1 % SAME locations used in second display
    xPos2 = xPos1;
    yPos2 = yPos1;
    cxPos2 = cxPos1;
    cyPos2 = cyPos1;
elseif condition == 2
%     secondDispCenters = datasample(nCenters(~ismember(nCenters,firstDispCenters)),nStim,'Replace',false); % randomly pick stimuli locs within quadrant
%     if iscolumn(secondDispCenters)
%         secondDispCenters = secondDispCenters';
%     end
%     
%     for l = 1:setSize/2
%         xPos2(l) = xGrid(secondDispCenters(l));
%         yPos2(l) = yGrid(secondDispCenters(l));
%     end 
%     for k = 1:(nStim - setSize/2)
%         cxPos2(k) = xGrid(secondDispCenters(k+setSize/2));
%         cyPos2(k) = yGrid(secondDispCenters(k+setSize/2));
%     end
    if setSize == 2
        pos = randi(3,1);
        xPos2 = cxPos1(pos);
        yPos2 = cyPos1(pos);
        cxPos2 = cxPos1;
        cxPos2(pos) = xPos1;
        cyPos2 = cyPos1;
        cyPos2(pos) = yPos1;
    else
        xPos2 = cxPos1;
        yPos2 = cyPos1;
        cxPos2 = xPos1;
        cyPos2 = yPos1;    
    end
else
%     firstDispCenters = datasample(nCenters,nStim,'Replace',false); % randomly pick stimuli locs
%     if iscolumn(firstDispCenters)
%         firstDispCenters = firstDispCenters'; % transpose if it's a col bc it doesn't work with nCenter = 1
%     end
    for l = 1:setSize
        xPos1(l) = xGrid(firstDispCenters(l));
        yPos1(l) = yGrid(firstDispCenters(l));
    end
    if setSize < nStim
        for k = 1:(nStim - setSize)
            cxPos1(k) = xGrid(firstDispCenters(k+setSize)); 
            cyPos1(k) = yGrid(firstDispCenters(k+setSize));
        end
    else
        cxPos1 = NaN(1);
        cyPos1 = NaN(1);
    end
    xPos2 = NaN(1);
    yPos2 = NaN(1);
    cxPos2 = NaN(1);
    cyPos2 = NaN(1);
end

jitter = round(prefs.stimSize*.25);

rAng = datasample(0:360,length(1:nStim),'Replace',true);
rAngT1 = rAng(1:length(xPos1));
rAngP1 = rAng(length(xPos1)+1:end);

if condition == 1
    rAngT2 = rAngT1;
    rAngP2 = rAngP1;
elseif condition == 2
    if setSize == 2
        rAngT2 = rAngP1(pos);
        rAngP2 = rAngP1;
        rAngP2(pos) = rAngT1;
    elseif setSize == 4
        rAngT2 = rAngP1;
        rAngP2 = rAngT1;
    end
end

xPos1 = round(xPos1 + (jitter*cosd(rAngT1)));
yPos1 = round(yPos1 + (jitter*sind(rAngT1)));

if condition ~= 3
    xPos2 = round(xPos2 + (jitter*cosd(rAngT2)));
    yPos2 = round(yPos2 + (jitter*sind(rAngT2)));
    
    cxPos1 = round(cxPos1 + (jitter*cosd(rAngP1)));
    cyPos1 = round(cyPos1 + (jitter*sind(rAngP1))); 
    
    cxPos2 = round(cxPos2 + (jitter*cosd(rAngP2)));
    cyPos2 = round(cyPos2 + (jitter*sind(rAngP2)));    
end

% for l = 1:length(xPos1)
%     if condition == 1
%         rAng = datasample(0:360,length(xPos1),'Replace',true); 
%         xPos1(l) = round(xPos1(l) + (jitter*cosd(rAng(l))));
%         yPos1(l) = round(yPos1(l) + (jitter*sind(rAng(l))));
% 
%         xPos2(l) = round(xPos2(l) + (jitter*cosd(rAng(l))));
%         yPos2(l) = round(yPos2(l) + (jitter*sind(rAng(l))));
%     elseif condition == 2
%         rAng = datasample(0:360,length(xPos1),'Replace',true);
%         xPos1(l) = round(xPos1(l) + (jitter*cosd(rAng(l))));
%         yPos1(l) = round(yPos1(l) + (jitter*sind(rAng(l))));
% 
%         rAng = datasample(0:360,length(xPos1),'Replace',true);
%         xPos2(l) = round(xPos2(l) + (jitter*cosd(rAng(l))));
%         yPos2(l) = round(yPos2(l) + (jitter*sind(rAng(l))));
%     else
%         rAng = datasample(0:360,length(xPos1),'Replace',true); 
%         xPos1(l) = round(xPos1(l) + (jitter*cosd(rAng(l))));
%         yPos1(l) = round(yPos1(l) + (jitter*sind(rAng(l))));
%     end
% end
% for l = 1:length(cxPos1)
%     if condition == 1
%        rAng = datasample(0:360,length(cxPos1),'Replace',true); 
%        cxPos1(l) = round(cxPos1(l) + (jitter*cosd(rAng(l))));
%        cyPos1(l) = round(cyPos1(l) + (jitter*sind(rAng(l)))); 
%     
%        cxPos2(l) = round(cxPos2(l) + (jitter*cosd(rAng(l))));
%        cyPos2(l) = round(cyPos2(l) + (jitter*sind(rAng(l))));
%     elseif condition == 2
%        rAng = datasample(0:360,length(cxPos1),'Replace',true);
%        cxPos1(l) = round(cxPos1(l) + (jitter*cosd(rAng(l))));
%        cyPos1(l) = round(cyPos1(l) + (jitter*sind(rAng(l))));
% 
%        rAng = datasample(0:360,length(cxPos1),'Replace',true);
%        cxPos2(l) = round(cxPos2(l) + (jitter*cosd(rAng(l))));
%        cyPos2(l) = round(cyPos2(l) + (jitter*sind(rAng(l))));
%     else
%        rAng = datasample(0:360,length(cxPos1),'Replace',true); 
%        cxPos1(l) = round(cxPos1(l) + (jitter*cosd(rAng(l))));
%        cyPos1(l) = round(cyPos1(l) + (jitter*sind(rAng(l))));
%     end
% end

end
%-------------------------------------------------------------------------
function pix = deg2pix(deg,vDist,pixSize)
% Convert degrees of visual angle to pixels for easy specification of
% stimulus size for PsychToolbox. Returns the size of a stimulus in
% pixels:
%
% INPUTS:
% deg: desired stim size in degrees of visual angle.
% vDist: viewing distance.
% pxSize: pixel size (in same units as viewing distance).
rad = deg2rad(deg); % convert visual angle from degrees to radians
sz = vDist*tan(rad); % size of stimulus (in same units as vDist and pxSize)
pix = round(sz/pixSize); % convert to pixels
end
%-------------------------------------------------------------------------
function win = openWindow(p) % open up the window! 

if numel(Screen('Screens')) > 1
    win.screenNumber = 1;
else
    win.screenNumber = 0;
end

%   p.windowed = 0; %%% 1  == smaller screen for debugging; 0 === full-sized screen for experiment
if p.windowed
    Screen('Preference', 'SkipSyncTests', 1);
    x_size=  1024; y_size = 768; %resolution, wxh
    [win.onScreen,rect] = Screen('OpenWindow', win.screenNumber, [128 128 128],[0 0 x_size y_size],[],[],[]);
    win.screenX = x_size;
    win.screenY = y_size;
    win.screenRect = [0 0 x_size y_size];
    win.centerX = (x_size)/2; % center of screen in X direction
    win.centerY = (y_size)/2; % center of screen in Y direction
    win.centerXL = floor(mean([0 win.centerX])); % center of left half of screen in X direction
    win.centerXR = floor(mean([win.centerX win.screenX])); % center of right half of screen in X direction
        % % Compute foreground and fixation rectangles
    win.foreRect = round(win.screenRect./1.5);
    win.foreRect = CenterRect(win.foreRect,win.screenRect);
    % compute eccentricity of color wheel
    win.colWheelRect = [0 0 350 350];
    win.colWheelRect = CenterRect(win.colWheelRect,win.screenRect);
else
%     Screen('Preference', 'SkipSyncTests', 1);
%     Screen('Preference','VisualDebugLevel',0);
    [win.onScreen,rect] = Screen('OpenWindow', win.screenNumber, [128 128 128],[],[],[],[]); %defaults to whole screen
    [win.screenX, win.screenY] = Screen('WindowSize', win.onScreen); % check resolution
    win.screenRect  = [0 0 win.screenX win.screenY]; % screen rect
    win.centerX = win.screenX * 0.5; % center of screen in X direction
    win.centerY = win.screenY * 0.5; % center of screen in Y direction
    win.centerXL = floor(mean([0 win.centerX])); % center of left half of screen in X direction
    win.centerXR = floor(mean([win.centerX win.screenX])); % center of right half of screen in X direction
    % % Compute foreground and fixation rectangles
%     win.foreRect = round(win.screenRect./1.5);
    win.foreRect = [0 0 967 720]; % 15 x 11.3 degrees of visual angle
    win.foreRect = [0 0 966 966]; % 15 x 11.3 degrees of visual angle
    win.foreRect = [0 0 766 766]; % 15 x 11.3 degrees of visual angle
    win.foreRect = CenterRect(win.foreRect,win.screenRect); % center the smaller, first rect in the second rect
    % compute eccentricity of color wheel
    win.colWheelRect = [0 0 350 350];
    win.colWheelRect = CenterRect(win.colWheelRect,win.screenRect);

    HideCursor; % hide the cursor since we're not debugging
end

% Screen('BlendFunction', win.onScreen, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA', [1 1 1 0]); % for transparency?

% basic drawing and screen variables
win.black    = BlackIndex(win.onScreen);
win.white    = WhiteIndex(win.onScreen);
win.gray     = mean([win.black win.white]);
win.darkgray = (win.gray+win.black)/2;
% win.gray     = win.darkgray;

win.foreColor = win.gray;

win.fontsize = 24;

win.colors_9 = [255, 0, 0;
0, 255, 0;
0, 0, 255;
255, 255, 0;
0, 255, 255];

win.cGrey = [191,191,191];

% win.colors_9 = [255 0 0; ... % red
%     0 255 0; ...% green
%     0 0 255; ...% blue
%     255 255 0; ... % yellow
%     255 0 255; ... % magenta
%     0 255 255; ... % cyan
%     255 255 255; ... % white
%     1 1 1; ... %black
%     255 128 0]; % orange!

% Set up a rect for the stimtrak probe dot to appear in!
win.probeRect = [win.screenX-40 10 win.screenX-10 40]; %%%% Upper right Corner!!!!!

% make a dummy call to GetSecs to load the .dll before we need it
dummy = GetSecs; clear dummy;
end
%-------------------------------------------------------------------------
%  CHANGE PREFERENCES!
%-------------------------------------------------------------------------
function prefs = getPreferences
% Design conditions
prefs.numBlocks = 20;
prefs.nTrialsPerCondition = 3;
prefs.conditions = [1 2 3]; % 1 = Same, 2 = Different, 3 = Simultaneous
prefs.setSizes = [2 4]; % 2,4
prefs.probeDisplay = [1 2]; % 1 = first display, 2 = second display
prefs.change = [0 1]; %0 = same color, 1 = different color

% Timing
prefs.refreshCycle =  1/120; % refresh rate in hertz
prefs.interStimInterval = .350;
prefs.retentionInterval = .900;
prefs.stimulusDuration = .150;
prefs.breakLength = .1*5; % number of minutes for block

% Stimulus size & positions
prefs.vDist = 75; % viewing distance (cm)
prefs.px = 0.0277; % pixel size (cm) 1920 X 1080
prefs.stimSize = deg2pix(2,prefs.vDist,prefs.px);
prefs.fixationSize = deg2pix(0.10,prefs.vDist,prefs.px);
prefs.circleSize = sqrt((prefs.stimSize^2)/pi);

% Randomize trial order of full factorial design order
prefs.fullFactorialDesign = fullfact([length(prefs.conditions), ...
    length(prefs.setSizes),...
    length(prefs.probeDisplay),...
    length(prefs.change),...
    prefs.nTrialsPerCondition]);

% Total number of trials in each fully-crossed block.
prefs.numTrials = size(prefs.fullFactorialDesign,1);

end