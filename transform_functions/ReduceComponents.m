function componentMask = ReduceComponents(binaryMask, componentMode)
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

    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    global pltCount;

    % Reduce size since our component counting implementation is quite
    % expensive
    componentMask = logical(resize(binaryMask, size(binaryMask)./2));

    %get biggest and center component
    [~, biggest, ~, center, centerCorners] = CountComponents(componentMask);

    if strcmp(componentMode, 'biggest')
        componentMask = biggest;
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(componentMask); title('Biggest Component');
        end
    elseif strcmp(componentMode, 'center')
        componentMask = center;
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            hold on;
            imshow(componentMask); title('Most Centered Component');
            plot(polyshape(centerCorners), 'EdgeColor', 'green', 'FaceColor', 'green');
            hold off;
        end
    end

    %resize to previous size
    componentMask = logical(resize(componentMask, size(binaryMask)));
end