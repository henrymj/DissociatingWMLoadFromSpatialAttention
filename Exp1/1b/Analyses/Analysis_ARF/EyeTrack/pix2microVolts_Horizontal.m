%function to take in each time point of the eye tracker gaze data, and to
%convert the X Y coordinates first to degrees of visual angle and then to
%microvolts so we can plot the data on the same scale as the VEOG and HEOG
%for comparrison. 


function [degH,microV_H] = pix2microVolts_Horizontal(pixH,vDist) %pix is num pixels, vDist is dist from screen
%position - fixation?
%Height and Width of Benq in CM

vDist = 72; % we'll start with our average distance from the screen for now 
%ideally we want this to update trial by trial or moment by moment. 

width = 53.5; 
height = 30;
%Resolution of BENQ
pixelsH = 1920; 
pixelsV = 1080;
%Pixel Size In CM
pixSizeH = width/pixelsH;
pixSizeV = height/pixelsV;

pixHfromFix = pixH-(pixelsH/2);

%Degrees for Horizontal
rad = atan((pixHfromFix)/vDist*pixSizeH);
degH = rad2deg(rad);
%Convert to Microvolts = deg * 20 microvolts
microV_H = degH*16;

end

