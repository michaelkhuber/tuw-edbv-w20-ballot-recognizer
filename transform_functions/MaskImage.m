function [maskedImage, pltCount] = MaskImage(img, pltCount)
% MASKIMAGE  is used to compute the edges between a distinct foreground
% (e.g. white ballot paper) and its background. 
% It does this by computing the image's gradient edges and then
% reducing the amount of components, only leaving the biggest components in the mask.
% Three steps are performed:
%   - The gradient magnitudes of the input image are computed for each
%   pixel. Then, a threshold is applied to the image, leaving only pixels
%   that have a higher magnitude than mean(magnitudes) + 1.0 *
%   std(magnitudes) where magnitudes are the gradient magnitudes of all
%   pixels from the image. Thus, only the most significant 32% edges
%   remain.
%   - Performs a Component Reduction on the resulting gradient mask with
%   strength 1.0 (see function ReduceComonents for more info).
%   - Dilates the resulting mask to thicken the edges.
%
% Author:
%   Richard Binder
%
% Source:
%   Self
%
% Inputs:
%   img:              The image to mask
%
% Output:
%   maskedImage:    resulting masked image

        global showPlot;
        global savePlot;
        global pltM;
        global pltN;
       

		% Create a gradient magnitude mask
		[gradMag, ~] = imgradient(img);
		gradThreshold = mean(gradMag(:)) + 1.0 * std(gradMag(:));
		gradMask = (gradMag > gradThreshold);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(gradMask); title('Grad Mask');
        end
        
		% Find all the connected compoments & remove small ones
        componentMask = ReduceComponents(gradMask, 1.0);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(componentMask); title('Component Reduction');
        end

        DilationMask = componentMask;
%         se = strel('octagon',12);
%         DilationMask = imdilate(DilationMask, se);
%         
%         if(showPlot || savePlot) 
%             subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
%             imshow(DilationMask); title('Dilation Mask');
%         end
        
        maskedImage = DilationMask;
        
end
