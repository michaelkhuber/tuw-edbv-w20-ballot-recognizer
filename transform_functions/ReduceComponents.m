function componentMask = ReduceComponents(binaryMask, strength)
% REDUCECOMPONENTS reduces the amount of components based on their size of a binary image mask. It
% removes all components that have size less than mean(AllSizes) + strength
% * std(AllSizes) where AllSizes is the list of all component sizes.
%
% Author:
%   Richard Binder
%
% Inputs:
%   binaryMask:                a binary image mask
%   strength:                   a value multiplied with standard deviation, e.g. if
% 1.0, then the smallest 68% components are removed
%
% Output:
%   componentMask:    resulting binary image mask with reduced amount of components

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
        
%         ccSizes = zeros(cc.NumObjects,1);
%         for i = 1 : cc.NumObjects
%             currCC = cc.PixelIdxList{i};
%             ccSizes(i) = size(currCC, 1);
%         end
%         
% 		ccSizeThreshold = mean(ccSizes(:)) + strength*std(ccSizes(:));
%         maxComponent = [];
%         
%         for i = 1 : cc.NumObjects
%             currCC = cc.PixelIdxList{i};
%             if size(currCC, 1) > size(maxComponent, 1)
%                 maxComponent = currCC;
%             end
%             if size(currCC, 1) < ccSizeThreshold
%                 componentMask(currCC) = 0;
%             end
%         end
%         
%         componentMask(maxComponent) = 1;
end