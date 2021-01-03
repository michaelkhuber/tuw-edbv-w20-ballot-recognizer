function [circle_centers, circle_radii] = FindCircles(im, radii)
% FINDCIRCLES Performs circular hough tranform on image, finds peaks and returns found circles
%
% Author: 
%   Michael Huber
%
% Sources:
%   Conceptually based on Pedersen, S. J. K. (2007). Circular hough transform. 
%   Aalborg University, Vision, Graphics, and Interactive Systems, 123(6).
%   
%   Implementation is based on David Young (2020). Hough transform for circles 
%   (https://www.mathworks.com/matlabcentral/fileexchange/26978-hough-transform-for-circles), 
%   MATLAB Central File Exchange. Retrieved December 27, 2020.
%
% Inputs:
%   im:             Image to find circles in
%   radii:          Radii of circles to find in image
% Outputs:
%   circle_centers: Nx2 matrix with each entry being x and y coordinate of circle found
%   circle_radii:   N vector with each entry being radius of circle found
%
% Uses:
%   CircularHough, CircularHoughPeaks

    % Histogram Equalization (maybe not neccessary)
    % im = adapthisteq(im);
    
    % Find edges in image
    im = edge(im, 'canny');
    
    % Create radii range
    radii_range = radii(2):-1:radii(1);

    % Compute Hough transform
    h = CircularHough(im, radii_range);
    
    % Find peaks in Hough space and output circles
    [circle_x, circle_y, circle_radii] = CircularHoughPeaks(h, radii_range);
    circle_centers = [circle_x, circle_y];
end
    
function h = CircularHough(im, radii)
% CIRCULARHOUGH Performs ciruclar hough transform for circles with given radii on binary image im
%   and returns accumulator array h.
%
% Author:
%   Michael Huber
%
% Source:
%   Implementation is based on David Young (2020). Hough transform for circles 
%   (https://www.mathworks.com/matlabcentral/fileexchange/26978-hough-transform-for-circles), 
%   MATLAB Central File Exchange. Retrieved December 27, 2020.
%
% Inputs:
%   im:     Binary image to find circles in
%   radii:  Radii of circles to find
%
% Output:
%   h:      Accumulator array containing votes in hough space
%
% Uses: RasteredCircle

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
% CIRCULARHOUGHPEAKS Finds the peaks in accumulator array h and returns the peak values as circle
%   center coordinates and their corresponding radii.
%
% Author:
%   Michael Huber
%
% Source:
%   Implementation is based on David Young (2020). Hough transform for circles 
%   (https://www.mathworks.com/matlabcentral/fileexchange/26978-hough-transform-for-circles), 
%   MATLAB Central File Exchange. Retrieved December 27, 2020.
%
% Inputs:
%   h:          Accumulator array obtained by Circular Hough Transform
%   radii:      Radii of circles to find
%
% Output:
%   x_out:      x-coordinates of circle centers
%   y_out:      y-coordinates of circle centers
%   radii_out:  radii of circles

    SENSITIVITY = 0.4;
    
    % Find local maxima in accumulator
    h_max = imregionalmax(h);
    h_temp = h_max;
    
    % Increase sensitivity until low number of maxima is found
    do = true;
    while do
        % Define SENSITIVITY% of maximum value as threshold
        threshold = SENSITIVITY * max(h, [], 'all');
        
        % Remove maxima below threshold
        h_temp = h_max & (h>=threshold);
        
        % Define number of circles found
        num_circles_found = length(find(h_temp));
        
        % Break or increase sensitivity if too many circles found
        if num_circles_found < 100
           break; 
        else
           SENSITIVITY = SENSITIVITY + 0.05;
        end
    end
    h_max = h_temp;

    % Get peaks
    peak_indices = find(h_max);
    [y, x, radii_indices] = ind2sub(size(h), peak_indices);
    peaks = [x'; y'; radii(radii_indices)];

    % Output results
    x_out = peaks(1,:)';
    y_out = peaks(2,:)';
    radii_out = peaks(3,:)';
end
    
function [x, y] = RasteredCircle(r)
% RASTEREDCIRCLE Generate rasterized circle using simple circle rastering algorithm based on 
%   rounding circle coordinates
%
% Author:
%   Michael Huber
%
% Source:
%   Implementation is based on David Young (2020). Hough transform for circles 
%   (https://www.mathworks.com/matlabcentral/fileexchange/26978-hough-transform-for-circles), 
%   MATLAB Central File Exchange. Retrieved December 27, 2020.
%
% Inputs:
%   r:  Radius of circle to generate
%
% Output:
%   x:  x-coordinates of circle
%   y:  y-coordinates of circle

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
