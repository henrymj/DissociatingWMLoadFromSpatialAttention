function [cloud_idx, protected_idx] = gen_cloud_idx(cloud, other_cloud_idx, other_protected_idx, stimulus)
    
    % get possible indices to sample from
    nBins = length(cloud.bins);
    possible_idx = NaN(nBins*stimulus.nDotsPerBin, 1);
    protected_first_idx = NaN(nBins, 1);
    protected_last_idx = NaN(nBins, 1);
    for bini = 1:nBins
        bin_start_idx = stimulus.binStarts(cloud.bins(bini));
        curr_idx =  bin_start_idx:bin_start_idx+stimulus.nDotsPerBin-1;
        possible_idx(1+(bini-1)*stimulus.nDotsPerBin:bini*stimulus.nDotsPerBin) = curr_idx;
        
        % note an idx in the first and last slice of each bin
        % to reserve for the second cloud
        protected_first_idx(bini) = randsample(curr_idx(1:stimulus.nRows), 1);
        protected_last_idx(bini) = randsample(curr_idx(end-stimulus.nRows+1:end), 1);
    end
    
    protected_idx = [protected_first_idx, protected_last_idx];
%     protected_idx = reshape([protected_first_idx, protected_last_idx], [], 1);
    cloud_idx = NaN(cloud.nDots, 1);
    
    if isempty(other_cloud_idx)
        % First cloud:
        % Must use a dot in the first and last slice of the
        % bins it spans. Those dots must not come from the protected set.
        % Remaining dots are sampled uniformly, avoiding the protected set.
        cloud_idx(1) = randsample(setdiff(possible_idx(1:stimulus.nRows), protected_idx), 1);
        cloud_idx(2) = randsample(setdiff(possible_idx(end-stimulus.nRows+1:end), protected_idx), 1);
        cloud_idx(3:end) = randsample(setdiff(possible_idx, protected_idx), cloud.nDots-2);
    else
        % Second cloud:
        % For first and last slice, check for a reserved dot.
        % Use if it exists.
        % If if does not exist, sample from the first and last slice
        % indices, while avoiding those used by the first cloud.
        % for the remaining dots, sample uniformly while avoiding those
        % used by the first dot cloud.
        
        % check for a reserved dot in the first slice, use if it exists
        if ~isempty(intersect(possible_idx(1:stimulus.nRows), other_protected_idx))
            intersection = intersect(possible_idx(1:stimulus.nRows), other_protected_idx);
            assert(length(intersection)==1);
            cloud_idx(1) = intersection;
        else
            cloud_idx(1) = randsample(setdiff(possible_idx(1:stimulus.nRows), other_cloud_idx), 1);
        end
        
        % check for a reserved dot in the last slice, use if it exists
        if ~isempty(intersect(possible_idx(end-stimulus.nRows+1:end), other_protected_idx))
            intersection = intersect(possible_idx(end-stimulus.nRows+1:end), other_protected_idx);
            assert(length(intersection)==1);
            cloud_idx(2) = intersection;
        else
            cloud_idx(2) = randsample(setdiff(possible_idx(end-stimulus.nRows+1:end), other_cloud_idx), 1);
        end
        cloud_idx(3:end) = randsample(setdiff(possible_idx, other_cloud_idx), cloud.nDots-2);
        
    end
end