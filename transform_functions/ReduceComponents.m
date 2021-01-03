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

        componentMask = binaryMask;
		cc = bwconncomp(binaryMask, 4);
        
        maxComponent = [];
        for i = 1 : cc.NumObjects
            currCC = cc.PixelIdxList{i};
            if size(currCC, 1) > size(maxComponent, 1)
                maxComponent = currCC;
            end
        end
        
        componentMask(:) = 0;
        componentMask(maxComponent) = 1;
end