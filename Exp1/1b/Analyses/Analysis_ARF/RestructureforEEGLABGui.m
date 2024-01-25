%Script to restructure EEG data from(Trials,chans,timepoints) to (Chans,Timepoints,Trials) for plotting with the EEG Lab plot function
%This script also concatenates the EEG data and Eye tracking data (if we have it for
%the subject) along with the logical indices from our auto rejection:
%EEG: Noise and Blocking          EyeTracking: Parsers Blinks and Saccades
%

%The script saves 2 files:
%1. A restructured .mat file (subj#_restruct_for_arf).
%2. An EEG lab file for manual marking (subj#MarkingComplete.set and .fdt
clear all

subjects = [38];
dbstop if error

    % add the eyetracking root to the path!!1
    addpath([pwd,'/EyeTrack'])

for s = 1:length(subjects)


sn = subjects(s);
%add eyetrack channels or don't 0 = no eye tracking file 1 = eye tracking
%file
eyetrack= 0;
eyeTrackPlotChans = 6; % 2 for eyetracker and 4 to plot fixation lines
noEyeTrackFixChans = 2; % add 2 "fixation lines" for reference with EOG
root = [pwd,'/EegData/',char(num2str(sn)),'/',];
%name some files
%Data restructure Name (what we're generating here)
rName = [root,num2str(sn),'_restruct_for_arf.mat'];
%pull in the compiled EEG data file
dName = [root,num2str(sn),'_EEG.mat'];
%pull in the eye tracking file
eName = [root,num2str(sn),'_eye_seg.mat'];
load(dName);
if eyetrack ==1
    nchans = size(erp.arfDat.baselined,2)+eyeTrackPlotChans; %size EEG plus 3 for eyetracker and 3 to plot fixation lines
    load(eName); %load the segmented eye track file
    nTrialEye = size(eyeData.epoched,1);
    nTimePtsEye = size(eyeData.epoched,3);
    nTrialEEG = size(erp.arfDat.baselined,1);

    %%% just to make it easier, downsample eyeData here!
    eyeData.epoched = eyeData.epoched(:,:,1:2:end);
else
    nchans = size(erp.arfDat.baselined,2)+noEyeTrackFixChans;
    nTrialEye = 0; nTrialEEG = 0;
end% drop out stim track from plotting (after checking it out)
%we only want to do the restructure step once.
if nTrialEye ~=nTrialEEG
    sprintf('We dont see eye to eye. Check  your EEG and EYE segment trial numbers')
else
    if exist(rName)
        sprintf('Already Did This Step. It is time to move on!')
    else
        %preallocate a matrix for the restructured data
        if eyetrack == 1

            tmpeeg = erp.arfDat.baselined(:,:,:);
            restructured = nan(size(tmpeeg,2)+eyeTrackPlotChans,size(tmpeeg,3),size(tmpeeg,1));
        else
            tmpeeg = erp.arfDat.baselined(:,:,:);
            restructured = nan(size(tmpeeg,2)+noEyeTrackFixChans,size(tmpeeg,3),size(tmpeeg,1));
        end

        for ch =1:size(tmpeeg,2)
            for tr = 1:size(restructured,3)
                restructured(ch,:,tr) = squeeze(tmpeeg(tr,ch,:));
            end
        end
        %Eye track data is structured differently
        if eyetrack ==1
            % Allocate variables for shifting the channels so they appear on
            % the butterfly plot
            HEOGShift = 100;
            HEyeShift = 200;
            VEOGShift = -200;
            VEyeShift = -100;
            StimTrackShift = 400;
            size(eyeData.epoched,1);

            for tr = 1:size(restructured,3)
                %Convert pixels to microvolts for Horiz Eye Movements
                [degV,microV] = pix2microVolts_Horizontal(squeeze(eyeData.epoched(tr,1,:)));
                restructured(nchans-5,:,tr) = microV;
                %plot horizontal lines so we can compare actual traces to a
                %horizontal
                %fixation reference We have to plot these as their own channels
                %in EEG lab. It's clunky but it works.
                restructured(nchans-3,:,tr) = HEOGShift;
                restructured(nchans-2,:,tr) = HEyeShift;
            end

            for tr = 1:size(restructured,3)
                %Convert pixels to microvolts for Vert Eye Movements
                [degV,microV] = pix2microVolts_Vertical(squeeze(eyeData.epoched(tr,2,:)));
                restructured(nchans-4,:,tr) = microV;
                %plot horizontal lines so we can compare actual traces to a
                %vertical
                %fixation reference  We have to plot these as their own channels
                %in EEG lab. It's clunky but it works.
                restructured(nchans-1,:,tr) = VEOGShift;
                restructured(nchans,:,tr) = VEyeShift;
            end
            %shift heog so we can see it on the butterfly plot
            restructured(nchans-8,:,:) = restructured(nchans-8,:,:)+HEOGShift;
            %shift veog
            restructured(nchans-7,:,:) = restructured(nchans-7,:,:)+VEOGShift;
            %shift stim track to top
            restructured(nchans-6,:,:) = restructured(nchans-6,:,:)+StimTrackShift;
            %shift horizontal Eye Track
            restructured(nchans-5,:,:) = restructured(nchans-5,:,:)+HEyeShift;
            %shift vertical Eye Track
            restructured(nchans-4,:,:) = restructured(nchans-4,:,:)+VEyeShift;
        else
            %Does the same thing just doesn't add the eye tracker or eye
            %tracker fixation chans.
            HEOGShift = 100;
            VEOGShift = 200;
            StimTrackShift = 400;
            for tr = 1:size(restructured,3)
                restructured(nchans-1,:,tr) = HEOGShift;
                restructured(nchans,:,tr) = VEOGShift;
            end
            restructured(nchans-4,:,:) = restructured(nchans-4,:,:)+HEOGShift;
            restructured(nchans-3,:,:) = restructured(nchans-3,:,:)+VEOGShift;
            restructured(nchans-2,:,:) = restructured(nchans-2,:,:)+ StimTrackShift;
        end
        save(rName,'restructured')
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

        % Here are the things you'll have to change! pnts = time points in your
        % file,  x min is the pre trial interval so the values look sensible!
        EEG = pop_importdata('dataformat','matlab','nbchan',nchans,'data',rName,'srate',500,'pnts',976,'xmin',-.400);
        EEG = eeg_checkset( EEG );
        %import auto reject
        %load the auto reject blocking variable into the rejconst variable so we
        %get a marked color indicator for blockingCURRENTSET
        EEG.reject(:).rejjp = (erp.arf.noise | erp.arf.dropout);
        %%% Add an adidtional one for blinks.....

        EEG.reject(:).rejconst = erp.arf.drift; %%% honestly just picked out a random ass name with "reject" in it sinc I didn't know the naming scema... KA
        %%% Add an adidtional one for blinks.....

        %load the auto reject noise variable into the rejthresh variable so we
        %get a marked color indicator for blocking
%         EEG.reject(:).rejkurt = erp.arf.noise;

        %load the PARSER BLINK function into eegLab Variable rejection structure
        if eyetrack ==1

            blinkSum = eyeData.calculatedBadEyeVals + erp.arf.blink;
            eMoveSum = eyeData.calculatedEMove + erp.arf.eMove;


            blinkSum(blinkSum>1) = 1; % make all 0's and 1's again.
            eMoveSum(eMoveSum>1) = 1; % make all 0's and 1's again.


            EEG.reject(:).rejthresh = blinkSum;

            %load the auto reject noise variable into the rejthresh variable so we
            %get a marked color indicator for PARSER BLINKS
            EEG.reject(:).rejfreq = erp.arf.eMove; %%% eog reject
            EEG.reject(:).rejkurt = eyeData.calculatedEMove;   %%%% parser reject

        else

            EEG.reject(:).rejthresh = erp.arf.blink;

            EEG.reject(:).rejfreq = erp.arf.eMove;
        end

        EEG = eeg_checkset( EEG );
        ename = [num2str(sn),'MarkingComplete_postTrial'];
        EEG = pop_saveset( EEG, 'filename',ename,'filepath',root);
    end
end
%Make sure to manually save dataset as the specified name after you finish
%the arf. Otherwise you'll lose the work you did.

end
