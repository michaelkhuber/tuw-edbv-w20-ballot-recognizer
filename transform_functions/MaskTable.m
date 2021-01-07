function maskedImage = MaskTable(img)
% MASKIMAGE2  is used to compute the edges of the ballot paper
% table. Only the edges of the table itself will remain, all others (edges
% between background and ballot paper, edges of ballot title "Amtlicher
% Stimmzettel") are removed.
% It performs the same three steps as the function MASKIMAGE,
% but changes some threshold values and adds two more steps, making it five as a whole:
%   - A threshold is applied to the image, leaving only pixels that have a higher gradient 
%   magnitude than max(magnitudes) * 0.07 where magnitudes are the gradient 
%   magnitudes of all pixels from the image. This is a form of edge
%   detection.
%   - A Component Reduction is performed on the resulting gradient mask 
%   (see function ReduceComponents for more info). This ensures
%   that only the table edges remain, as they make up the biggest component.
%   - The Convex Hull of the resulting image mask (the table) is computed
%   - The edges of that hull (by using gradient again) are computed. This
%   leaves us with only the outter boundaries of the table.
%   - The resulting mask is dilated to thicken the edges for more reliable 
%   Hough Transformation.
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
    gradThreshold = max(gradMag(:)) * 0.08;
    gradMask = (gradMag > gradThreshold);

    if(showPlot || savePlot) 
        pltCount = pltCount + 1; subplot(pltM, pltN, pltCount);
        t = 'Grad Mask';
        imshow(gradMask); title([num2str(pltCount), '. ', t]);
    end

    % Find all the connected compoments & get biggest one
    componentMask = ReduceComponents(gradMask, 'center');

    ConvexHull = bwconvhull(componentMask);
    ConvexHull(:,size(ConvexHull,2)) = 0.0;
    ConvexHull(:,1) = 0.0;
    ConvexHull(size(ConvexHull,1),:) = 0.0;
    ConvexHull(1, :) = 0.0;

    if(showPlot || savePlot) 
        pltCount = pltCount + 1; subplot(pltM, pltN, pltCount);
        t = 'Convex Hull';
        imshow(ConvexHull); title([num2str(pltCount), '. ', t]);
    end

    [gradMag, ~] = imgradient(ConvexHull);
    ConvexHull = gradMag > 0.00001;

    structure = se('octagon',6);
    ConvexHull = dilate(ConvexHull, structure);

    if(showPlot || savePlot) 
        pltCount = pltCount + 1; subplot(pltM, pltN, pltCount);
        t = 'Border Convex Hull';
        imshow(ConvexHull); title([num2str(pltCount), '. ', t]);
    end

    maskedImage = ConvexHull;
end