function [x, y, sizes, dot_colors, cloud1, cloud2] = gen_dot_info(stimulus, cloud1, cloud2)
    n_slices = stimulus.nSlices; % slices of the circle
    n_rows = stimulus.nRows; % rows of the grid
    min_dist  = stimulus.minEccentricity_pix; % start of closest row, in pixels
    max_dist = stimulus.maxEccentricity_pix;
    row_thickness = round((max_dist - min_dist)/n_rows);
    size_min = stimulus.minDotSize_pix;
    size_max = stimulus.maxDotSize_pix;
    slice_thickness = stimulus.sliceThickness;
  
    % dot size(s)
    sizes = repmat((size_min:(size_max-size_min)/(n_rows-1):size_max)', 1, n_slices);

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

    % cloud1
    [cloud1_idx, protected_idx] = gen_cloud_idx(cloud1, [], [], stimulus);
    cloud1.idx = cloud1_idx;
    
    
    % dot colors
    dot_colors = repmat(stimulus.placeholderColor, n_rows*n_slices, 1);
    dot_colors(cloud1.idx, :) = repmat(stimulus.dotColors.(cloud1.color), cloud1.nDots, 1);
    
    
    if isa(cloud2, 'struct')
        [cloud2_idx, ~] = gen_cloud_idx(cloud2, cloud1_idx, protected_idx, stimulus);
        cloud2.idx = cloud2_idx;
        dot_colors(cloud2.idx, :) = repmat(stimulus.dotColors.(cloud2.color), cloud2.nDots, 1);
    end

    % reshape to vectors for diplay, output
    x = reshape(x, [1, n_slices*n_rows]);
    y = reshape(y, [1, n_slices*n_rows]); 
    sizes = reshape(sizes, [1, n_slices*n_rows]);    
    dot_colors = dot_colors';
end