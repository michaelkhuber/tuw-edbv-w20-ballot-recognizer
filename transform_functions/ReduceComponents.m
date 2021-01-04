function componentMask = ReduceComponents(binaryMask)
% REDUCECOMPONENTS reduces the amount of components based on their size in a binary image mask. It
% removes all components such that only the biggest one is left behind.
%
% Author:
%   Richard Binder
%
% Inputs:
%   binaryMask:                a binary image mask
%
% Output:
%   componentMask:    resulting binary image mask with biggest component

        % Reduce size since our component counting implementation is quite
        % expensive
        componentMask = logical(resize(binaryMask, size(binaryMask)./2));
        
        %get biggest component
        [num_components, biggest, biggest_size] = CountComponents(componentMask);
        
        %resize to previous size
        componentMask = logical(resize(biggest, size(binaryMask)));
end