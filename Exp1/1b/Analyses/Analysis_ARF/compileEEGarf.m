function compileEEGarf

subs = [38];
root = pwd;

for i = 1:length(subs)

    dataRoot = ['/EegData/',num2str(subs(i)),'/'];
    %Load in the EEG Data so EEG and beh variables are saved together
    dRoot = [root,dataRoot];
    cd(dRoot)
    fName = [root,dataRoot,num2str(subs(i)),'_EEG.mat'];
    load(fName)

    %load the EEG Lab ARF file
    EEG = pop_loadset('filename',[num2str(subs(i)),'MarkingComplete_postTrial.set'],'filepath',dRoot);
    EEG = eeg_checkset( EEG );

    erp.arf.m_blockingNoiseDropout = EEG.reject(:).rejjp;
    erp.arf.m_drift = EEG.reject(:).rejconst;
    erp.arf.m_blink = EEG.reject(:).rejthresh;
    erp.arf.m_saccade_heog = EEG.reject(:).rejfreq;
    erp.arf.m_saccade_eyetrack = EEG.reject(:).rejkurt;
    erp.arf.manual = EEG.reject(:).rejmanual;

    % in case we forgot to save the file... (no manual !)
    if isempty(erp.arf.manual)
        fprintf('No Manual detected!! Remember to save dataset!!!')
    end
    for t =  1:size(erp.trial.data,1)
        erp.arf.artifactIndCleaned(t) = erp.arf.m_blockingNoiseDropout(t) | erp.arf.m_drift(t) | erp.arf.m_blink(t) | erp.arf.m_saccade_heog(t) | erp.arf.m_saccade_eyetrack(t) | erp.arf.manual(t);
    end
    erp.arf.proportion_arfs = sum(erp.arf.artifactIndCleaned) ./ length(erp.arf.artifactIndCleaned);
    erp.arf.trials_remaining = length(erp.arf.artifactIndCleaned) - sum(erp.arf.artifactIndCleaned);
    fprintf('\nProportion Arfs: %.2f \n',erp.arf.proportion_arfs);
    fprintf('Trials Remaining: %d \n',erp.arf.trials_remaining);

    %%% experiment-specific: check for the minimum number of trials in any
    %%% given conidtion
    numCond1 = sum(erp.trial.codes'==11 & erp.arf.artifactIndCleaned==0);
    numCond2 = sum(erp.trial.codes'==12 & erp.arf.artifactIndCleaned==0);
    numCond3 = sum(erp.trial.codes'==13 & erp.arf.artifactIndCleaned==0);
    numCond4 = sum(erp.trial.codes'==14 & erp.arf.artifactIndCleaned==0);
    numCond5 = sum(erp.trial.codes'==15 & erp.arf.artifactIndCleaned==0);
    numCond6 = sum(erp.trial.codes'==19 & erp.arf.artifactIndCleaned==0);
    numCond7 = sum(erp.trial.codes'==20 & erp.arf.artifactIndCleaned==0);
    numCond8 = sum(erp.trial.codes'==21 & erp.arf.artifactIndCleaned==0);
    numCond9 = sum(erp.trial.codes'==22 & erp.arf.artifactIndCleaned==0);
    numCond10 = sum(erp.trial.codes'==23 & erp.arf.artifactIndCleaned==0);


    erp.arf.num_per_cond = [numCond1,numCond2,numCond3,numCond4,numCond5,numCond6,numCond7,numCond8,numCond9,numCond10];
    erp.arf.min_num_per_cond = min(erp.arf.num_per_cond)';

        fprintf('NumCond1: %d \n',numCond1);
        fprintf('NumCond2: %d \n',numCond2);
        fprintf('NumCond3: %d \n',numCond3);
        fprintf('NumCond4: %d \n',numCond4);
        fprintf('NumCond5: %d \n',numCond5);
        fprintf('NumCond6: %d \n',numCond6);
        fprintf('NumCond7: %d \n',numCond7);
        fprintf('NumCond8: %d \n',numCond8);
        fprintf('NumCond9: %d \n',numCond9);
        fprintf('NumCond10: %d \n',numCond10);
        fprintf('Min Per Cond: %d \n',erp.arf.min_num_per_cond);

    save(fName,'erp','-v7.3');
    clear erp

    fprintf('\n File saved! \n');

end

cd(root)
