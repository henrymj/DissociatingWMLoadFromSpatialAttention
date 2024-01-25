%function to take in each time point of the eye tracker gaze data, and to
%convert the X Y coordinates first to degrees of visual angle and then to
%microvolts so we can plot the data on the same scale as the VEOG and HEOG
%for comparrison. 
% 
% written by Kirsten Adam  to calculate the change in degrees of visual
% angle that occurs over the course of a window length (calculated with our
% typical step function.... 
% http://www.en.eyebrainpedia.com/#!practical-considerations-vs-theoretical/c11s7

function deg = pix2deg_changeInPos(prePeakX,postPeakX,prePeakY,postPeakY,vDist) %pix is num pixels, vDist is dist from screen
%position - fixation?
%Height and Width of Benq in CM

vDist = 100; % we'll start with our average distance from the screen for now 
% ideally we want this to update trial by trial or moment by moment. 

%Dimensions of BENQ monitors
width = 53.5; 
height = 30;
%Resolution of BENQ
pixelsH = 1920; 
pixelsV = 1080;
%Pixel Size In CM
pixSizeH = width/pixelsH;
pixSizeV = height/pixelsV;


%% Starting position (degrees from fixation)
pixHfromFix = prePeakX-(pixelsH/2);
pixVfromFix = prePeakY-(pixelsV/2);
%Degrees for Horizontal
rad = atan((pixHfromFix/2)/vDist*pixSizeH)*2;
degH1 = rad2deg(rad);

%Degrees for Vertical
rad = atan((pixVfromFix/2)/vDist*pixSizeV)*2;
degV1 = rad2deg(rad);
%% Ending position (degrees from fixation)
pixHfromFix = postPeakX-(pixelsH/2);
pixVfromFix = postPeakY-(pixelsV/2);
%Degrees for Horizontal
rad = atan((pixHfromFix/2)/vDist*pixSizeH)*2;
degH2 = rad2deg(rad);

%Degrees for Vertical
rad = atan((pixVfromFix/2)/vDist*pixSizeV)*2;
degV2 = rad2deg(rad);

%% Distance between starting and ending poitn (degrees)

deg = sqrt((degH2-degH1).^2 + (degV2-degV1).^2);

end

