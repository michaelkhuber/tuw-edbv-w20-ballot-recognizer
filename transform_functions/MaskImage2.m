function maskedImage = MaskImage2(img)
% MASKIMAGE2  is used to compute the edges of the ballot paper
% table. Only the edges of the table itself will remain, all others (edges
% between background and ballot paper, edges of ballot title "Amtlicher
% Stimmzettel") are removed.
% It performs the same three steps as the function MASKIMAGE,
% but changes some threshold values and adds two more steps, making it five as a whole:
%   - The gradient magnitudes of the input image are computed for each
%   pixel. Then, a threshold is applied to the image, leaving only pixels
%   that have a higher magnitude than mean(magnitudes) - 0.2 *
%   std(magnitudes) where magnitudes are the gradient magnitudes of all
%   pixels from the image. Thus, approximately the most significant 60% edges
%   remain.
%   - Performs a Component Reduction on the resulting gradient mask with
%   strength 3.0 (see function ReduceComonents for more info). This ensures
%   that only the table edges remain, as they make up the biggest component.
%   - Dilates the resulting mask to thicken the edges.
%   - Computes the Convex Hull of the resulting image mask (the table)
%   - Computes the edges of that hull (by using gradient again). This
%   leaves us with only the outter edges of the table (as opposed to inner ones).
%   - Again dilates the resulting mask to thicken the edges.
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
    global pltCount;

    % Create a gradient magnitude mask
    [gradMag, ~] = imgradient(img);
    %gradThreshold = mean(gradMag(:)) - 0.2 * std(gradMag(:));
    gradThreshold = max(gradMag(:)) * 0.07;
    gradMask = (gradMag > gradThreshold);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
        imshow(gradMask); title('Grad Mask');
    end

    % Find all the connected compoments & remove small ones
    componentMask = ReduceComponents(gradMask);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
        imshow(componentMask); title('Component Reduction');
    end  

    ConvexHull = bwconvhull(componentMask);
    ConvexHull(:,size(ConvexHull,2)) = 0.0;
    ConvexHull(:,1) = 0.0;
    ConvexHull(size(ConvexHull,1),:) = 0.0;
    ConvexHull(1, :) = 0.0;

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        imshow(ConvexHull); title('Convex Hull');
    end

    [gradMag, ~] = imgradient(ConvexHull);
    ConvexHull = gradMag > 0.00001;

    structure = se('octagon',6);
    ConvexHull = dilate(ConvexHull, structure);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        imshow(ConvexHull); title('Border Convex Hull');
    end

    maskedImage = ConvexHull;
end