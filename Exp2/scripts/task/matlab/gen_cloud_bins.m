function bins = gen_cloud_bins(stimulus)
    cloud_nBins = randsample(stimulus.binSizes, 1);
    firstBin = randi([1, stimulus.nBins]);
    bins = firstBin + (1:cloud_nBins) - 1;
    % bin wraparound
    bins= mod(bins, stimulus.nBins);
    bins(bins==0) = stimulus.nBins;
end