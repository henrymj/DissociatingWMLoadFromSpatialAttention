function [x, y, sizes, colors] = gen_dot_info(cloud_1_nDots, cloud_2_nDots)
    n_bins = 8; % CTF bins
    n_slices = 64; % slices of the circle
    n_rows = 6; % rows of the grid
    row_thickness = 60; % arbitrary radial distance in pixels of each row of the grid
    min_dist  = 160; % start of closest row, in pixels
    size_min = 16;
    size_max = 19;
    bin_sizes = [1, 2]; % possible number of bins a cloud can occupy

    min_grey_per_bin = 6;
    max_grey_per_bin_wColors = 8;
    max_grey_per_bin_woColors = 48;

%     cloud_1_nDots = 10;
%     cloud_2_nDots = 20;

    color_map = containers.Map('KeyType','char','ValueType','any');
    color_names = {'red'; 'green'; 'blue'; 'yellow'; 'purple'; 'teal'; 'orange'};
    color_vals = [[255, 0, 0]; [0, 255, 0]; [0, 0, 255]; [255, 255, 0]; [255, 0, 255]; [0, 255, 255]; [255, 128, 0]];
    for colori = 1:length(color_names)
        color_map(char(color_names(colori))) = color_vals(colori, :);
    end

    cloud_1_color = char(randsample(color_names, 1));
    cloud_2_color = char(randsample(setdiff(color_names, cloud_1_color), 1));


    % constants from inputs
    slice_thickness = 2*pi / n_slices; % radians per slice
    nDots_per_bin = n_slices*n_rows/n_bins;
    bin_starts = 1:nDots_per_bin:nDots_per_bin*(n_bins-1) + 1;

    % dot size(s)
    sizes = repmat((size_min:(size_max-size_min)/(n_rows-1):size_max)', 1, 64);
%     sizes = size_min + rand(n_rows, n_slices)*(size_max-size_min);
    % pixel_sizes = sizes./10; %TODO: refine for displays

    % radial distances (used to constrain angles)
    distance_jitter = rand(n_rows, n_slices).*(row_thickness-sizes);
    distances = repmat(min_dist + (0:row_thickness:(n_rows-1)*row_thickness)', 1, n_slices) + distance_jitter;

    % angles
    initial_angles = 0:slice_thickness:2*pi-slice_thickness;
    angle_jitter = (rand(n_rows, n_slices)-.5).*(slice_thickness - (sizes./distances));
    angles = repmat(initial_angles, n_rows, 1) + angle_jitter;

    % final locations

    x = cos(angles).*distances;
    y = sin(angles).*distances;

    % clouds
    [cloud_1_idx, cloud_1_bins] = gen_cloud_idx(cloud_1_nDots, [], bin_sizes, n_bins, nDots_per_bin, bin_starts);
    [cloud_2_idx, cloud_2_bins] = gen_cloud_idx(cloud_2_nDots, cloud_1_idx, bin_sizes, n_bins, nDots_per_bin, bin_starts);

    % random control dots
    grey_indices = [];
    for bini = 1:n_bins
        possible_indices = bin_starts(bini):bin_starts(bini)+nDots_per_bin-1;
        if ismember(bini, cloud_1_bins) || ismember(bini, cloud_2_bins)
            n_grey = randi([min_grey_per_bin, max_grey_per_bin_wColors]);
            possible_indices = setdiff(setdiff(possible_indices, cloud_1_bins), cloud_2_bins);
        else
            n_grey = randi([min_grey_per_bin, max_grey_per_bin_woColors]);
        end
        grey_indices = cat(1, grey_indices, randsample(possible_indices, n_grey)');
    end

    % colors
    colors = repmat([166, 166, 166], n_rows*n_slices, 1);
    colors(cloud_1_idx, :) = repmat(color_map(cloud_1_color), cloud_1_nDots, 1);
    colors(cloud_2_idx, :) = repmat(color_map(cloud_2_color), cloud_2_nDots, 1);


    % reshape to vectors for diplay, output
    x = reshape(x, [1, n_slices*n_rows]);
    y = reshape(y, [1, n_slices*n_rows]); 
    sizes = reshape(sizes, [1, n_slices*n_rows]);
    
    % outputs
%     final_idx = cat(1, grey_indices, cloud_1_idx, cloud_2_idx);
%     x = x(final_idx);
%     y = y(final_idx);
%     sizes = sizes(final_idx);
%     colors = colors(final_idx, :)';
    colors = colors';
end