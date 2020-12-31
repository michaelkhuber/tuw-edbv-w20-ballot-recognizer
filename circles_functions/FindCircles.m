function [circle_centers, circle_radii] = FindCircles(im, radii)
    %%% Michael Huber
    
    % Convert to binary
    im = edge(im, 'canny');
    
    % Create radii range
    radii_range = radii(1):1:radii(2);

    % Compute Hough transform
    h = CircularHough(im, radii_range);
    
    % Find peaks in Hough space and output circles
    [circle_x, circle_y, circle_radii] = CircularHoughPeaks(h, radii_range);
    circle_centers = [circle_x, circle_y];
end
    
function h = CircularHough(im, radii)
% Get indices of non-zero pixels in binary image
    [im_y, im_x] = find(im);
    
    % Initialize accumulator array with padding to avoid out of bounds errors
    [size_y, size_x] = size(im);
    padding = ceil(max(radii));
    h = zeros(size_y+2*padding, size_x+2*padding, length(radii));
    h_size = size(h);
    
    % Generate rastered circles from radii
    circles_x = []; circles_y = []; circles_radii = [];
    for i = 1:length(radii)
       [circle_x, circle_y] = RasteredCircle(radii(i));
       circles_x = [circles_x circle_x];
       circles_y = [circles_y circle_y];
       circles_radii = [circles_radii repmat(i,1,length(circle_x))];
    end
    
    % Loop over image pixels
    for i = 1:length(im_x)
        % Place circles on pixel for votes
        indices1 = padding+im_y(i)+circles_y;
        indices2 = padding+im_x(i)+circles_x;
        indices3 = circles_radii;
    
        % Calculate linear indices of votes
        votes = sub2ind(h_size, indices1, indices2, indices3);
        
        % Vote in the accumulator
        h(votes) = h(votes)+1;
    end
    
    % Remove padding from accumulator array
    h = h(1+padding:end-padding, 1+padding:end-padding, :);
end

function [x_out, y_out, radii_out] = CircularHoughPeaks(h, radii)

    % Define 50% of maximum value as threshold
    threshold = 0.5 * max(h, [], 'all');

    % Find local maxima in accumulator
    h_max = imregionalmax(h);

    % Remove maxima below threshold
    h_max = h_max & (h>=threshold);

    % Get peaks
    peak_indices = find(h_max);
    [y, x, radii_indices] = ind2sub(size(h), peak_indices);
    peaks = [x'; y'; radii(radii_indices)];

    % Sort by strength if more than n results
    %if n < size(peaks,2)
    %    [~, indices_sorted] = sort(h(h_max), 'descend');
    %    indices_sorted = indices_sorted(1:n);
    %    peaks = peaks(:, indices_sorted);
    %end

    % Output results
    x_out = peaks(1,:)';
    y_out = peaks(2,:)';
    radii_out = peaks(3,:)';
end
    
function [x, y] = RasteredCircle(r)
    % Generate rasterized circle using simple circle rastering algorithm
    % based on rounding the values.
    
    % Get the x-coordinate of the point 45° in the upper right (NE) quadrant,
    % which therefore is the end of the first octant (NNE). This point is the 
    % one where x=y, therefore x = r/sqrt(2).
    x_end = round(r/sqrt(2));
    
    % Generate coordinates for the NNE octant
    x_oct = 0:x_end;
    y_oct = round(sqrt(r.^2-x_oct.^2));
    
    % Assemble NE quadrant by adding the mirrored octant
    x_quad = [x_oct y_oct(end:-1:2)]; 
    y_quad = [y_oct x_oct(end:-1:2)];
    % Assemble E circle half by adding the mirrored quadrant
    x_half = [x_quad y_quad];
    y_half = [y_quad -x_quad];
    % Assemble full circle by adding the mirrored half
    x = [x_half -x_half];
    y = [y_half -y_half];
    end
