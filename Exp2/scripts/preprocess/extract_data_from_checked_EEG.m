subs = {'007', '008'};
numsubs= length(subs);
sep = filesep;
baseDir = fileparts(fileparts(pwd));
destination = [baseDir, sep, 'DATA', sep, 'Preprocessed',sep];
summary_log = fopen('extraction_summary.txt','w');
eeglab

for isub = 1:numsubs
    checked_file = [baseDir, sep, 'DATA', sep, 'EEG', sep, subs{isub}, sep, [subs{isub}, '_DCD_checked.set']];
    EEG = pop_loadset(checked_file);
    
    % Summary
    fprintf(summary_log, ['\nRunning ', EEG.setname, '\n\n']);
    trial_counts = rej_summary(EEG);
    for x = 1:size(trial_counts, 1)
        fprintf(summary_log,'Condition %1.f: %1.f\n',trial_counts(x, 1),trial_counts(x, 2));
    end
    
    % Titles
    title = regexprep(EEG.setname, ' +', '_');
    xdata_filename = [destination, title, '_xdata.mat'];
    ydata_filename = [destination, title, '_ydata.mat'];
    idx_filename = [destination, title, '_artifact_idx.mat'];
    behavior_filename = [destination, title, '_behavior.csv'];
    info_filename = [destination, title, '_info.mat'];
    
    % Remove unwanted channels and save xdata
    num_chans = EEG.nbchan;
    all_chans = strings(num_chans,1);
    for chan = 1:num_chans
        all_chans(chan,:) = EEG.chanlocs(chan).labels;
    end
    chan_idx = ismember(all_chans,{'L_GAZE_X','L_GAZE_Y','R_GAZE_X','R_GAZE_Y','StimTrak','HEOG','VEOG','TP9','GAZE_X','GAZE_Y','L-GAZE-X','L-GAZE-Y','R-GAZE-X','R-GAZE-Y','GAZE-X','GAZE-Y'});

    xdata = EEG.data(~chan_idx,:,:);
    save(xdata_filename, 'xdata');
    
    
    % ydata
    num_trials = size(xdata,3);
    ydata = zeros(num_trials,1);
    for x=1:num_trials
        sorted_labels = sort(EEG.epoch(x).eventbinlabel);
        char_labels = char(sorted_labels(end));
        ydata(x,:) = str2double(char_labels(6:end-1));
    end
    
    save(ydata_filename, 'ydata');
    
    
    % Saving artifact index for indexing behavior file
    num_rows = size(EEG.event,2);
    all_trials = zeros(num_rows,1);
    for x = 1:num_rows
        all_trials(:,x) = EEG.event(x).bepoch;
    end
    checked_trials = unique(all_trials);
    
    unchecked_file = regexprep(checked_file, 'checked', 'unchecked');
    unchecked_EEG = pop_loadset(unchecked_file);
    unchecked_trials = (1:unchecked_EEG.trials)';
    artifact_idx = ismember(unchecked_trials,checked_trials);
    save(idx_filename,'artifact_idx')
    
    % Save copy of behavior csv
    behavior_file_old = [baseDir, sep, 'DATA', sep, 'Behavior', sep, 'p-',subs{isub}, '_task-changeDetection_Exp1.csv'];
    copyfile(behavior_file_old,behavior_filename);
    
    
    % Gather other info variables
    chan_labels = {EEG.chanlocs.labels}';
    chan_labels = char(chan_labels(~chan_idx));
    chan_x = [EEG.chanlocs.X];
    chan_y = [EEG.chanlocs.Y];
    chan_z = [EEG.chanlocs.Z];
    chan_x = chan_x(~chan_idx);
    chan_y = chan_y(~chan_idx);
    chan_z = chan_z(~chan_idx);
    sampling_rate = EEG.srate;
    times = EEG.times;
    
    load(regexprep(behavior_file_old, '.csv', '.mat'))
    unique_id = str2num(participant.ID);

    save(info_filename,'unique_id','chan_labels','chan_x','chan_y','chan_z','sampling_rate','times');

     
%     clear labels num_trials templabel x y checked_trials participant
    
end
% clear all