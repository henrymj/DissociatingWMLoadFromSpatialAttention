function erp = preprocessSettings(erp)
% This function specifies the settings that the preprocessing functions call on.

fprintf('loading preprocessing settings... \n')

erp.droppedElectrodes = {'FT9','FT10'};  % specificy unwanted electrodes (these are our 2 references)

erp.codes = [11:14,19:22];

% Timing for artifact rejection (times should be ABSOLUTE VALUES) 
erp.arfPreTime = 400; % msec prior to timelock
erp.arfPostTime = 1550; % msec post timelock

% Timing stuff for building waveforms (times should be absolute values) 
erp.preTime = erp.arfPreTime+500; % msecs prior to timelock built in extra 500ms for time freq analyses
erp.postTime = erp.arfPostTime+500; % msecs post timelock

% Window for baselining (time should be NEGATIVE IF PRE-=TIMELOCKING
erp.baseStart = -400;  %%% if using the whole time period, just use -erp.preTime
erp.baseEnd = 0;

% Noise threshold for artifact rejection
erp.arf.noiseThr = 150; % microvolts

% Threshold for drift
erp.arf.driftThr = 100; %microvolts

% Step function settings for channel drop out (main cap channels sometimes
% have step functions when they suddenly drop out! 
% do a wide window length to avoid catching alpha!! 
erp.arf.dropoutWin = 250; %ms
erp.arf.dropoutStep = 20; % ms
erp.arf.dropoutThr = 100; % microvolts

% Step function settings for blink rejection
erp.arf.blinkWin = 150; % ms
erp.arf.blinkStep = 10; % ms
erp.arf.blinkThr = 30; % microvolts
% Step function settings for horizontal eye movements 
erp.arf.eMoveWin = 150; % ms
erp.arf.eMoveStep = 10; % ms
erp.arf.eMoveThr = 20; %microvolts
