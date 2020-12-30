function [centers_out,radii_out] = SortCircles(centers_in,radii_in)
%SORTCIRCLES Summary of this function goes here
%   Detailed explanation goes here

% Pixel threshold, circle centers within 10px of each other are treated
% as equal
THRESHOLD = 20;

% Include first circle
centers_out = [centers_in(1,:)];
radii_out = [radii_in(1)];

% Loop through other circles
for i = 2:length(centers_in)
    % Check if circle is equal to already covered circle
    difference = abs(centers_out(:,2) - centers_in(i,2));
    not_in_group = difference > THRESHOLD;
    newgroup = all(not_in_group);
    
    % Include if no equal circle found
    if newgroup
        centers_out = [centers_out; centers_in(i,:)];
        radii_out = [radii_out; radii_in(i)];
    end
end


end

