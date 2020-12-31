function [centers_out,radii_out] = FilterCircles(centers_in,radii_in)
% FILTERCIRCLES Filters circles by removing all clustered circles that are similar to an already existing circle
%
% Author: 
%   Michael Huber
%
% Inputs:
%   centers_in:     Circle centers obtained by FindCircles
%   radii_in:       Radii of circles obtained by FindCircles
% Outputs:
%   centers_out:    Filtered circle centers
%   radii_out:      Filtered cicle radii

    global Y_THRESHOLD;
    global X_THRESHOLD;

    % Y threshold, circle centers within 20px of each other are treated
    % as equal
    Y_THRESHOLD = 20;

    % X-coordinate median and X threshold
    x_median = median(centers_in(:,1));
    X_THRESHOLD = 50;

    % Include first circle that is within median threshold
    start_index = 1;
    for i = 1:length(centers_in)
        if abs(centers_in(i,1)-x_median) < X_THRESHOLD
            centers_out = centers_in(i,:);
            radii_out = radii_in(i);
            start_index = i;
            break
        end
    end

    % Loop through other circles
    for i = start_index:length(centers_in)
        % Check if circle is equal to already covered circle
        difference = abs(centers_out(:,2) - centers_in(i,2));
        in_group = difference < Y_THRESHOLD;
        newgroup = all(not(in_group));
        
        % Include if no equal circle found
        if newgroup
            % Only include if circle is near x median
            if abs(centers_in(i,1)-x_median) < X_THRESHOLD
                centers_out = [centers_out; centers_in(i,:)];
                radii_out = [radii_out; radii_in(i)];
            end
        else
            % Use circle if radius is bigger
            matching_groups = find(in_group);
            if length(matching_groups) == 1
                if radii_out(matching_groups) < radii_in(i)
                    radii_out(matching_groups) = radii_in(i);
                    centers_out(matching_groups) = centers_in(i);
                end
                
            % Error if more than one matching group is found
            else
                error("More than one group matching circle found");
            end 
        end
    end


end

