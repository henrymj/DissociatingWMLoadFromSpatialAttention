function [degV,microV_V] = pix2microVolts_Vertical(pixV,vDist) %pix is num pixels, vDist is dist from screen
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

pixVfromFix = pixV-(pixelsV/2);
%Degrees for Vertical
rad = atan((pixVfromFix)/vDist*pixSizeV);
degV = rad2deg(rad);
%Convert to Microvolts = deg * 20 microvolts
microV_V = degV*16;

end

