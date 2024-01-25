function bins = get_bins_from_loc_width(loc, width, nBins)
    bins = loc + (1:width) - 1;
    % bin wraparound
    bins= mod(bins, nBins);
    bins(bins==0) = nBins;
end