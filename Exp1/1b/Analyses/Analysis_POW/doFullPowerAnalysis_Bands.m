function doFullPowerAnalysis_Bands

% Calculate phase angle across all channels, frequencies, and times
% try, matlabpool local, end

dbstop if error

subjects = [1,2,3,4,7,8,9,10,11,12,17,18,20,21,23,33,34,36,37,38];

nSubs = length(subjects);
root = pwd;

% load a sample file
fName = [root,'/EegData/',char(num2str(subjects(1))),'/',char(num2str(subjects(1))),'_EEG.mat'];
load(fName);

Fs = 500; %Sampling frequency

p.freqs = [8 12];  % theta, alpha, beta
p.freqNames = {'Alpha'};
p.time = erp.trial.times;
nFreqBands = size(p.freqs,1);
nPoints =  length(p.time);

for i = 1:nSubs
    % get the subject number
    subdata = char(num2str(subjects(i)));
    sn = subjects(i);
    fprintf('Subject:\t%d\n',sn)

    p.channels = erp.chanLabels;
    whichChans = ~ismember(p.channels,{'HEOG','VEOG','StimTrak'});
    p.channels = p.channels(whichChans);

    % load the subject's file
    root = pwd;
    fName = [root,'/EegData/',char(num2str(subjects(i))),'/',char(num2str(subjects(i))),'_EEG.mat'];
    load(fName)

    % figure out the number of trials and conditions
    nTrials = ceil(size(erp.trial.data,1));
    nChans = length(p.channels);

    tCells = nChans*nFreqBands*nTrials;
    cCnt = 1;

    % set up space to save
    cPower = NaN(nChans,nFreqBands,nTrials,nPoints);
        % loop through channels (except eye channels)
        for ch = 1:nChans
            data = erp.trial.data(:,whichChans,:);
            data = squeeze(data(:,ch,:));

            % show how far along we are in the analysis on printed out on
            % the screen
            if cCnt>1
                pc = round(cCnt/tCells,2)*100;
                fprintf('Percent Complete:\t%d\n',pc)
            end

            % loop through each frequency band
            for f = 1:nFreqBands

                low_cut_off = p.freqs(f,1); % cut off for the filter
                high_cut_off = p.freqs(f,2);

                % Filter Data
                % try the field trip bandpass function instead
                filtered_data = ft_preproc_bandpassfilter(data,Fs,[low_cut_off,high_cut_off]);

                nActualTrials = size(filtered_data,1);

                parfor t = 1:nActualTrials
                    [amp,inst_phase,cum_phase]= hilbert_amp_phase(filtered_data(t,:));
                    cPower(ch,f,t,:) = amp;

                    cCnt = cCnt+1;
                end;
            end;
        end;
    p.power = cPower;

    saveName = [root,'/EegData/',subdata,'/',subdata,'_Power_FreqBands_ft.mat'];
    save(saveName,'p','fName','-v7.3');

    clear cPower

end;
