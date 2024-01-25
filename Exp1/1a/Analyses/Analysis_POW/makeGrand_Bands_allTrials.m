
function makeGrand_Bands_allTrials
load allsubs_data.mat % load behavior

subjects = [1,2,4,5,6,7,8,9,11,12];
nSubs = length(subjects);

% conditions!!!
nTimes = 976;
nConds = 4;
nFreqs = 1;
condInd_ss = [2 2 4 4];
condInd_condition = [1 2 1 2]; % 1 = same, 2 = different
gpow.condLabels = {'SS2 Same','SS2 Different','SS4 Same','SS4 Different'};

% allocate some space to save stuff!!
gpow.frontal = NaN(nSubs,nFreqs,nConds,nTimes);
gpow.central = NaN(nSubs,nFreqs,nConds,nTimes);
gpow.posterior = NaN(nSubs,nFreqs,nConds,nTimes);

% load a random erp file so we have crap!!
    loadFile = [pwd,'/EEGData/',char(num2str(subjects(1))),'/',char(num2str(subjects(1))),'_EEG.mat'];
    load(loadFile);

    nChans = erp.nChans-3;
    gpow.allchans =  NaN(nSubs,nFreqs,nChans,nConds,nTimes);


% loop through subs and calculateeeee
for s = 1:length(subjects);

    loadFile = [pwd,'/EegData/',char(num2str(subjects(s))),'/',char(num2str(subjects(s))),'_Power_FreqBands_ft.mat'];
    load(loadFile);

    loadFile2 =  [pwd,'/EegData/',char(num2str(subjects(s))),'/',char(num2str(subjects(s))),'_EEG.mat'];
    load(loadFile2);

    fprintf('Subject %d loaded! \n',subjects(s))

    ssInd = allsubs.setSize(s,:);
    condInd = allsubs.condition(s,:);

    gpow.filtTimes = p.time;
    gpow.trialTimes = -erp.arfPreTime:2:erp.arfPostTime;
    gpow.chans = erp.chanLabels; chans = gpow.chans; gpow.chans = gpow.chans(1:30);

    fchans = {'Fp1','Fp2','Fz','F3','F4'}; fchans = ismember(chans,fchans);
    cchans = {'Cz','C3','C4','CP1','CP2'}; cchans = ismember(chans,cchans);
    pchans = {'O1','O2','Oz','PO7','PO8','PO3','PO4','P3','P4','Pz','P7','P8'}; pchans = ismember(chans,pchans);
    pchansL = {'O1','P7','P3','PO7','PO3'}; pchansL = ismember(chans,pchansL);
    pchansR = {'O2','P8','P4','PO8','PO4'};  pchansR = ismember(chans,pchansR);

    for f = 1:nFreqs
        % set up temporary data
        tempData = squeeze(p.power(:,f,:,:));
        % only takethe relevant times, not the extra 500 ms
        tempData = tempData(:,:,ismember(p.time,-400:1550));

        if length(erp.arf.artifactIndCleaned) > size(tempData,2)
            erp.arf.artifactIndCleaned = erp.arf.artifactIndCleaned(1:size(tempData,2));
        end
        % replace artifacts with NaN's!!!!
        tempData(:,erp.arf.artifactIndCleaned==1,:) = NaN; % hand cleaned is "artifactIndCleaned"

        %%% replace incorrect trials with NaN's.....(optional) don't do for now
        %%% b/c wastes a lot of trials!!!
        %     tempData(accInd==0,:,:) = NaN;

        for c = 1:nConds

            behaviorIndex = ssInd==condInd_ss(c) & condInd==condInd_condition(c);

            if length(behaviorIndex)>size(tempData,2)
                behaviorIndex = behaviorIndex(1:size(tempData,2));
            end

            % FRONTAL gpow
            gpow.frontal(s,f,c,:) = squeeze(nanmean(nanmean(tempData(fchans,behaviorIndex==1,:),1),2));

            % CENTRAL gpow
            gpow.central(s,f,c,:) = squeeze(nanmean(nanmean(tempData(cchans,behaviorIndex==1,:),1),2));

            % POSTERIOR gpow
            gpow.posterior(s,f,c,:) = squeeze(nanmean(nanmean(tempData(pchans,behaviorIndex==1,:),1),2));

            %%% all chans
            for ch = 1:nChans
                gpow.allchans(s,f,ch,c,:) = squeeze(nanmean(nanmean(tempData(ch,behaviorIndex==1,:),1),2));
            end

        end % end condition loop

    end % end freqs loop

    fprintf('Subject %d Complete \n',subjects(s))

end % end subject loop

fprintf('Grand file made!! Savinggggggggggg \n')
save('grand_power_filteredBands_alltrials.mat','gpow','p','subjects','-v7.3');
fprintf('Saved!! \n')
