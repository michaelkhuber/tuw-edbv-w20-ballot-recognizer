function [lines, lines2] = HoughTransform(maskedImage, peakThreshold)
%HOUGHTRANSFORM finds horizontal and vertical lines in an image mask by
%using the hough transformation for lines.
% 
% Author:
%   Richard Binder
%
% Source:
%   All the Hough functions that were implemented are based on the Matlab
%   functions "hough", "houghpeaks" and "houghlines". They take similar
%   inputs and give similar outputs.
%   Helpful for the exact implementation was prior knowledge about Hough
%   Transformation and the documentation of the given matlab functions on
%   mathworks.com.
%   Also the (very) short demo on Hough Transformation at
%   https://www.mathworks.com/matlabcentral/fileexchange/52183-simple-demo-of-hough-transform-implementation
%   was helpful.
%
% Input:
%   maskedImage:    the image mask in which to find lines
%
% Output:
%   lines:      vertical lines that were found in the image mask
%   lines2:    horizontal lines that were found in the image mask

    global showPlot;
    global savePlot;
    global pltCount;

    rhoPrecision = 10;

    maxRho = sqrt(size(maskedImage,1)^2 + size(maskedImage,2)^2);
    maxRho = ceil(maxRho);
    rhoInterval = [-maxRho : rhoPrecision : maxRho];

    % Use Hough transform to find vertical lines
    thetaInterval = [-45.0 : 0.25 : 44.75];
    [H, thetaInterval, rhoInterval] = Hough(maskedImage, rhoInterval, thetaInterval);
    P = Houghpeaks(H, peakThreshold);
    lines = Houghlines(maskedImage, thetaInterval, rhoInterval, P);

    % Use Hough transform to find horizontal lines
    thetaInterval2 = [-90.0 : 0.25 : -44.75];
    thetaInterval2 = [thetaInterval2, 45.0 : 0.25 : 89.75];
    [H2, thetaInterval2, rhoInterval2] = Hough(maskedImage, rhoInterval, thetaInterval2);
    P2 = Houghpeaks(H2, peakThreshold);
    lines2 = Houghlines(maskedImage, thetaInterval2, rhoInterval2, P2);
    
    [H, thetaInterval, lines, lines2, P, P2] = mergeHoughSpaces(H, H2, thetaInterval, thetaInterval2, P, P2, lines, lines2);

    if(showPlot || savePlot)
         plotHoughTransform(H, thetaInterval, rhoInterval, P, P2);
    end
end





function [ H, thetaInterval, rhoInterval ] = Hough( maskedImage, rhoInterval, thetaInterval )
%HOUGH creates an accumulator matrix from a masked image. Each row in the
% matrix refers to a rho value, and each column refers to a theta value,
% where rho and theta are parameters for a line in hough representation.
% Each element in the accumulator matrix is the number of pixels on a particular line.
%
% Author:
%   Richard Binder
%
% inputs:
%   maskedImage:    the masked image
%   rhoInterval:          the interval of possible rho values for lines in
%   the masked image.
%   thetaInterval:          the interval of possible theta values for lines in
%   the masked image.
%
% outputs:
%   H:                          the accumulator array
%   thetaInterval:        same as input thetaInterval
%   rhoInterval:           same as input rhoInterval

    H = zeros(length(rhoInterval),length(thetaInterval));
    [y, x] = find(maskedImage);

    for thetaIndex = 1:length(thetaInterval)
        theta = thetaInterval(thetaIndex);
        rho = x*cosd(theta) + y*sind(theta);
        rhoIndices = getRhoIndex(rho, rhoInterval);
        thetaIndices =repmat(thetaIndex, [length(rhoIndices) 1]);
        H = H + accumarray([rhoIndices, thetaIndices],1,size(H));
    end
end

function rhoIndex = getRhoIndex(rho, rhoInterval)
% Author:
%   Richard Binder
    maxRho = max(rhoInterval);
    minRho = min(rhoInterval);
    rhoIndex = (rho + abs(minRho) ) / ( abs(minRho) + maxRho ) * length(rhoInterval);
    rhoIndex = floor(rhoIndex) + 1;
end





function P = Houghpeaks(H, peakThreshold)
%HOUGHPEAKS finds local maxima in a matrix H, whilst neglacting those that
%are smaller than peakThreshold * max(H(:))
%
% Author:
%   Richard Binder
%
% inputs:
%   H:      a matrix in which to find peaks
%
% outputs:
%   P:      the peaks in H, given as indices

    peaks = ordfilt2(H, 9, ones(3, 3));
    peaksBinary = (H == peaks) & (H > 0);
    
    threshold = peakThreshold * max(H(:));
    peaksBinary = peaksBinary & (peaks > threshold);
    
    [rho, theta] = find(peaksBinary);
    P = [rho, theta];
end





function lines = Houghlines(maskedImage, thetaInterval, rhoInterval, P)
%HOUGHLINES returns lines within a masked image in hough representation 
% based on the intervals of possible rho and theta values and their given indices
%
% Author:
%   Richard Binder
%
% inputs:
%   maskedImage:      the masked image
%   thetaInterval:         interval of possible theta values
%   rhoIntevral:            interval of possible rho values
%   P:                           list of lines, given as indices of
%   thetaInterval and rhoInterval
%
% outputs:
%   lines:      a list of structures. each structure is a line with
%   properties "rho", "theta", "point1", "point2"

    lines = [];
    for i = [1 : size(P, 1)]
        peak = P(i,:);
        line.rho = rhoInterval(peak(1));
        line.theta = thetaInterval(peak(2));
        
        [line.point1, line.point2] = findEndpoints(maskedImage, line, rhoInterval);
        
        lines = [lines, line];
    end
end





function [point1, point2] = findEndpoints(maskedImage, line, rhoInterval)
%FINDENDPOINTS finds the first and last endpoint of a line in a masked image
%given in hough representation. 
%
% Author:
%   Richard Binder
%
% inputs:
%   maskedImage:      the masked image
%   line:                        the line of which to find endpoints
%   rhoInterval:            the possible rho values for lines
%
% outputs:
%   point1:         the first endpoint
%   point2:         the second endpoint

    [y, x] = find(maskedImage);
    
    rho = x*cosd(line.theta) + y*sind(line.theta);
    rhoIndices = getRhoIndex(rho, rhoInterval);
    lineRhoIndex = getRhoIndex(line.rho, rhoInterval);
    
    coords = [x, y];
    linePixels = coords(rhoIndices == lineRhoIndex, :);
    linePixels = sortrows(linePixels, [1 2]);
    
    point1 = linePixels(1,:);
    point2 = linePixels(length(linePixels),:);
end







function [H, thetaInterval, lines, lines2, P, P2] = mergeHoughSpaces(H, H2, thetaInterval, thetaInterval2, P, P2, lines, lines2)
%MERGEHOUGHSPACES merges two hough spaces with different thetaIntervals into one, assuming that the first
%hough space goes between the two halfs of the second hough space. It also
%fixes an incorrect seperation of horizontal and vertical lines by assuming
%that lines in one of those groups have a theta difference of less than 45
%degrees.
% 
% Author:
%   Richard Binder
%
% Source:
%   Self
%
% Input:
%   H:                              the first Hough space
%   H2:                            the second Hough space
%   thetaInterval:            the theta interval of the first Hough space
%   thetaInterval2:          the theta interval of the second Hough space
%   P:                              the peaks of the first Hough space
%   P2:                            the peaks of the second Hough space
%   lines:                         the lines of the first Hough space
%   lines2:                       the lines of the second Hough space
%
% Output:
%   H:                              the merged Hough space
%   thetaInterval:            the merged theta interval
%   lines:                         vertical lines
%   lines2:                       horizontal lines
%   P:                              peaks of vertical lines
%   P2:                            peaks of horizontal lines

    %merge H and H2
    firstHalf = 1 : (floor(size(H2,2)/2)+1);
    secondHalf = (floor(size(H2,2)/2)+2) : size(H2,2);
    
    P(:,2) = P(:,2) + firstHalf(end);
    P2(P2(:,2) > firstHalf(end), 2) = P2(P2(:,2) > firstHalf(end), 2) + size(H, 2);
    H = [ H2(:, firstHalf), H, H2(:, secondHalf) ];
    thetaInterval = [ thetaInterval2(:, firstHalf), thetaInterval, thetaInterval2(:, secondHalf) ];
    
    %fix potentially incorrect seperation of horizontal and vertical lines
    allLines = [lines lines2];
    allP = [P; P2];
    
    verticalTheta = lines(1).theta;
    lines = [];
    lines2 = [];
    P =[];
    P2 = [];
    for i = 1:length(allLines)
        l = allLines(i);
        p = allP(i,:);
        
        diff1 = mod(l.theta -  verticalTheta, 180.0);
        diff2 = mod(verticalTheta - l.theta, 180.0);
        
        if min(diff1, diff2) < 45.0
            lines = [lines l];
            P = [P; p];
        else
            lines2 = [lines2 l];
            P2 = [P2; p];
        end
    end
end







function plotHoughTransform(H, thetaInterval, rhoInterval, P, P2)
%PLOTHOUGHTRANSFORM plots the accumulator matrix of a hough 
% transformation and its peaks.
% The plot also shows peaks of vertical and horizontal lines.
%
% Author:
%   Richard Binder
%
% Source:
%   The code is based on the hough transform visualization demo on
%   mathworks.com:
%   https://www.mathworks.com/help/images/ref/hough.html
%
% inputs:
%   H:                          the accumulator matrix for horizontal lines
%   thetaInterval:        the interval of possible theta values for horizontal lines
%   rhoInterval:           the interval of possible rho values for horizontal lines
%   P:                          the peaks in H for vertical lines, given as indices of
%   thetaInterval and rhoInterval.
%   P2:                          the peaks in H for horizontal lines, given as indices of
%   thetaInterval and rhoInterval.

    global pltM;
    global pltN;
    global pltCount;
    
    %plot Vertical Hough transformation
    subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
    imshow(imadjust(rescale(H)),'XData',thetaInterval,'YData',rhoInterval, 'InitialMagnification','fit');
     axis on, axis normal, hold on;
     xlabel('\theta'), ylabel('\rho');
    plot(thetaInterval(P(:,2)),rhoInterval(P(:,1)),'s','color','green'); %plot peaks
    plot(thetaInterval(P2(:,2)),rhoInterval(P2(:,1)),'s','color','blue'); %plot peaks
    colormap(gca,hot);
     title('Hough Transform - Vertical and Horizontal Peaks');
    hold off;
end