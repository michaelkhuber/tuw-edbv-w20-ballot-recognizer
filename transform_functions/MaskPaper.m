function maskedImage = MaskPaper(img)
% MASKIMAGE  is used to compute the edges between a distinct foreground
% (e.g. white ballot paper) and its background. 
% It does this by computing the image's gradient edges and then
% reducing the amount of components.
% Three steps are performed:
%   - A threshold is applied to the image, leaving only pixels that have a higher gradient 
%   magnitude than max(magnitudes) * 0.2 where magnitudes are the gradient 
%   magnitudes of all pixels from the image. This is a form of edge
%   detection.
%   - A Component Reduction is performed on the resulting gradient mask 
%   (see function ReduceComponents for more info). This ensures
%   that only the paper edges remain, as they make up the biggest component.
%   - The resulting mask is eroded to reduce the thickness of the edges to
%   get a more reliable Hough Transformation.
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
    gradThreshold = max(gradMag(:)) * 0.2;
    gradMask = (gradMag > gradThreshold);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
        imshow(gradMask); title('Grad Mask');
    end

    % Find all the connected compoments & remove small ones
    componentMask = ReduceComponents(gradMask, 'biggest');

    structure = se('octagon',9);
    ErodeMask = erode(componentMask, structure);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
        imshow(ErodeMask); title('Erosion Mask');
    end

    maskedImage = ErodeMask;
end
