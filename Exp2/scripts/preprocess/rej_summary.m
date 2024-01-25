function [t] = rej_summary(EEG)
    
    ydata = zeros(EEG.trials,1);
    for x=1:EEG.trials
        sorted_labels = sort(EEG.epoch(x).eventbinlabel);
        char_labels = char(sorted_labels(end));
        ydata(x,:) = str2double(char_labels(6:end-1));
    end
    
    unique_conds = unique(ydata);
    nconds = length(unique_conds);
    
    t = zeros(nconds,2);
    for x=1:nconds
        t(x, :) = [unique_conds(x), sum(ydata==unique_conds(x))];
    end