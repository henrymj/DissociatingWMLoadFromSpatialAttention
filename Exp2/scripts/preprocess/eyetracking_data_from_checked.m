subs = {'06','07','12','13','14','16','17','20','21','22','23','26','27','28','29','30','31','33','34','35','36','37','39','40','41','42','43','44','45','46','47','48','49','50'};
% subs = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31'};
% subs = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
experiment = '1801';
numsubs= length(subs);
destination = ['..\analysis\data\',experiment,'\'];
summary_log = fopen('summary.txt','w');
% eeglab

for isub = 1:numsubs
        processed_file = ['..\data\',experiment,'\',subs{isub},'\wst',experiment,'_',subs{isub},'_checked.set'];
        EEG = pop_loadset(processed_file);
        
        %Summary
        fprintf(summary_log, ['\nRunning ', EEG.setname, '\n\n']);
        trial_counts = rej_summary(EEG);
        
        %if subject has fewer than 200 in any condition, reject
        if min(trial_counts) < 200
            continue
        end
        
        %Titles
        eyetracking_filename = [destination, title, '_eyetracking.mat'];
        
        % Remove unwanted channels and save xdata
        num_chans = EEG.nbchan;
        all_chans = strings(num_chans,1);
        for chan = 1:num_chans
            all_chans(chan,:) = EEG.chanlocs(chan).labels;
        end
        chan_idx = ismember(all_chans,{'L_GAZE_X','L_GAZE_Y','R_GAZE_X','R_GAZE_Y','HEOG','VEOG','GAZE_X','GAZE_Y','GAZE-X','GAZE-Y'});

        eyetracking_data = EEG.data(chan_idx,:,:);
        eyetracking_labels = {all_chans{chan_idx}};
%         save(eyetracking_filename, 'eyetracking_data', 'eyetracking_labels');
        
        clear labels num_trials templabel EEG eyetracking_data
end
"DATA EXTRACTION COMPLETE"