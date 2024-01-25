clear all;
close all;
Screen('CloseAll');
KbName('UnifyKeyNames');

% Stimulus Details
stimulus.fixationSize_DVA = .2;        % DVA = degrees of visual angle
stimulus.fixationColour = 0;
stimulus.regionHeight_DVA = 12;

stimulus.minEccentricity_DVA = 2.5;
stimulus.maxEccentricity_DVA = 8;
stimulus.minDotSize_DVA = .2;
stimulus.maxDotSize_DVA = .35;

% Equipment
equipment.viewDist_MM = 800;           % ! MEASURE !
equipment.pixPerMM = 3.6; %4.3;        %3.6      % ! MEASURE ! using MeasureDpi function

equipment.greyVal = .5;      
equipment.gammaVals = [1 1 1];

% Calculate Pixel Space 
equipment.MMperDeg = (equipment.viewDist_MM/2)*tan(deg2rad(2*stimulus.maxEccentricity_DVA))/stimulus.maxEccentricity_DVA;
equipment.PixPerDeg = equipment.pixPerMM*equipment.MMperDeg;


stimulus.minEccentricity_pix = round(stimulus.minEccentricity_DVA*equipment.PixPerDeg);
stimulus.maxEccentricity_pix = round(stimulus.maxEccentricity_DVA*equipment.PixPerDeg);

% Psychtoolbox Set-Up  
AssertOpenGL;

screenID = 1; %max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');

[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[xCenter, yCenter] = RectCenter(winRect);

displaying = true;
while displaying
    colorReply = GetEchoString(ptbWindow, 'colVal (r,g,b): ', xCenter-100, yCenter, 0);
    currColor = cellfun(@str2num, strsplit(colorReply, ','));
    Screen('Flip', ptbWindow);
    Screen('FillRect', ptbWindow, currColor/255, [xCenter-stimulus.maxEccentricity_pix, yCenter-stimulus.maxEccentricity_pix, xCenter+stimulus.maxEccentricity_pix, yCenter+stimulus.maxEccentricity_pix]);
    Screen('DrawDots', ptbWindow, [xCenter, yCenter]', 5, [0, 0, 0], [], 2);
    Screen('Flip', ptbWindow);
    clear KbWait;
    WaitSecs(4)
    [secs, keyCode, deltaSecs] = KbWait(-1);
    if keyCode(KbName('RightArrow'))==1
        displaying = false;
        break
    end
end

sca;