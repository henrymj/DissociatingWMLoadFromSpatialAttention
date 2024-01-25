

function [blinkInd,markDegMove] = calc_eyeMove_dist_eyeTrack(eyeData)

%% Extract epoched X & Y Coordinates

x = squeeze(eyeData.epoched(:,1,:)); % x coordinates (pixels)
y = squeeze(eyeData.epoched(:,2,:)); % y coordinates (pixels)
s = squeeze(eyeData.epoched(:,3,:)); % pupil size 
d = squeeze(eyeData.epoched(:,4,:)); % distance from monitor 

%% remove nonsensical values where the eye tracker lost the eye because it was closed!
% mark trials where eye was lost as blinks!

nT = size(x,1);

for t = 1:nT
    tempX = x(t,:);
    tempY = y(t,:);
    
    tempX(tempX>2000 | tempX<0) = NaN;
    tempY(tempY>1500 | tempY<0) = NaN;
    
    x(t,:) = tempX;
    y(t,:) = tempY;
end

%% Use the nonsensical values as a proxy for blinks or where the eye tracker lost the eyes!

blinkInd = sum(isnan(x),2) > 0 | sum(isnan(y),2) > 0;
blinkInd = blinkInd';

%% Calculate max X distance and max Y distance, figure out if it violates our maximum eye displacement
%% Use a step function just like we would with EOG data. This will help us avoid false alarming
%  to those annoying eyelash flutters, etc. where the tracker blips but
%  it's not really an eye movement!


screenX = 1920; % pixels 
screenY = 1080; % pixels

markDegMove = NaN(1,nT);
for t = 1:nT
    
    if eyeData.mode == 1 %%%% desktop mode
        viewDist = eyeData.DistanceToScreen; 
    else
        viewDist = nanmean(d(t,:)) + eyeData.DistanceToCamera; % cm
    end
    degMove = step_eyeTrack(x(t,:),y(t,:),eyeData.stepSize,eyeData.winSize,eyeData.maxDeg,eyeData.rateAcq,viewDist);
    
    markDegMove(t) = artDetect(degMove);
end

end 


