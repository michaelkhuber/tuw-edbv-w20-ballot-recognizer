function transformed = Transform(im, resultName, step)
% TRANSFORM finds a ballot paper in the given input image and
% transforms the image such that the ballot paper becomes a (near) perfect rectangle. 
% Then crops the image to only the ballot paper. If step == 1, then the image is
% transformed and cropped to the ballot paper, if step == 2, then the image
% is transformed and cropped to the ballot table.
%
% Author:
%   Richard Binder
%
% Source:
%   Self
%
% Inputs:
%   im:                        The input image
%   resultName:        The name of the resulting plot file if savePlot is true
%   step:                     Wether to transform to ballot paper (1) or table (2)
%
% Output:
%   transformed:    The transformed and cropped ballot paper (or table)

    % Global Variables used for plotting
    global showPlot;
    global savePlot;
    global ballotFilename;
    global pltM;
    global pltN;
    global pltCount;
    
    showPlot = true; %show the plots in a figure
    savePlot = false; %save the plots as an image in a subfolder (expensive operation)
    ballotFilename = resultName;
    pltM = 3;
    pltN = 5;
    pltCount = 1;

    if step ~= 1 && step ~= 2
        error("Invalid transform step");
    end

    %resize image to have 2000 pixels in height without deforming it
    newSize = [2000, ceil(size(im,2)/size(im,1) * 2000)];
    im = resize(im,newSize);

    if(showPlot)
        f = figure(1);
        clf('reset');
    elseif(savePlot)
        f = figure('visible','off','Renderer', 'opengl', 'Position', [10 10 2000 1000]);
    end

    if(showPlot || savePlot)
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        imshow(im); title('Original');
    end

    % prepare the image
    if step == 1
        preparedImage = Prepare(im, 20.0, true);
    elseif step == 2
        preparedImage = Prepare(im, 5.0, false);
    end

    % mask the image
    if step == 1
        maskedImage = MaskImage(preparedImage);
    elseif step == 2
        maskedImage = MaskImage2(preparedImage);
    end

    % get lines from hough transformation of the masked image
    [lines, lines2] = HoughTransform(maskedImage);

    % find intersections between vertical and horizontal lines
    intersections = findIntersections(im, lines, lines2);

    % find the convex hull from all the intersection points
    k = convhull(intersections);
    hull = [];
    hull(:,:) = intersections(k,:);

    % reduces the hull to a quadrangle
    quadrangle = findCorners(hull);

    % orders the corners as top-left, top-right, bottom-right, bottom-left
    [corners, middlePoints, middleImagePoints] = orderCorners(im, quadrangle);

    if(showPlot || savePlot) 
        plotLinesAndCorners(im, lines, lines2, intersections, hull, corners, middlePoints, middleImagePoints);
    end

    % Transforms the image such that the ballot corners make up a
    % (near) perfect rectangle.
    imNew = HomographyTransformation(im, corners);
    
    if step == 2
        % Correct Table rotation if it is off
        imNew = correctTableRotation(imNew);
        
        %resize ballot table to have its original template size
        newSize = [2150,3500];
        imNew = resize(imNew,newSize);
    end

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
        if step == 1
            t = 'Cropped';
        else
            t = 'Rotated, Resized & Cropped'; 
        end
        imshow(imNew); title(t);
    end

    % Save plots as an image
    if(savePlot)
        print(f,sprintf("resources/results/Step2_Transform%i/Transform%i_%s", step, step, ballotFilename),'-dpng','-r700'); 
        clf(f);
        clear f;
    end

    transformed = imNew;
end





function intersections = findIntersections(im, lines, lines2)
    %FINDINTERSECTIONS finds intersection points between lines and lines2. 
    % If an intersetion point is outside the borders defined by the image im, then that point is ignored.
    %
    % Author:
    %   Richard Binder
    %
    % Source:
    %   Self
    %
    % Inputs:
    %   im:           The input image
    %   liines:       Set of lines defined by hough parameters rho and theta
    %   liines2:     Second set of lines defined by hough parameters rho and theta
    %
    % Output:
    %   intersections:    The intersection points
    
    imH = size(im, 1);
    imW = size(im, 2);
    
    intersections = [];
    for i = 1:length(lines)
      for j = 1:length(lines2)
        r1 = lines(i).rho;  t1 = lines(i).theta;
        r2 = lines2(j).rho; t2 = lines2(j).theta;

        x = (r1*sind(t2)-r2*sind(t1))/(cosd(t1)*sind(t2)-sind(t1)*cosd(t2));
        y = (r1*cosd(t2)-r2*cosd(t1))/(sind(t1)*cosd(t2)-cosd(t1)*sind(t2));
        if (isnan(x) || isnan(y) || x < 1-200 || x > imW+200 || y < 1-200 || y > imH+200)
            continue;
        end
        intersections = [intersections; x, y];
      end
    end
end

function quadrangle = findCorners(hull)
    %FINDCORNERS finds 4 corners points of hull, assuming that hull is a polygon
    %that resembles a quadrangle with 4 corners.
    %
    % Author:
    %   Richard Binder
    %
    % Source:
    %   This function uses a (way simpler) variation of the Ramer–Douglas–Peucker -
    %   Algorithm. It iteratively removes a point from the hull that is
    %   closest to the line segment between its two neighbouring points 
    % (making it the most irrelevant point of the hull), until only 4 points are left.
    %
    % Inputs:
    %   hull:       a set of ordered points making up a convex hull
    %
    % Output:
    %   quadrangle:    a set of 4 corner points, making up a quadrangle that resembles the
    %   input hull
    
    quadrangle = hull;
    while (length(quadrangle) > 4)
        minDist = inf;
        n = length(quadrangle);
        removeIndex = 0;

        for i=1:n
            pt = [quadrangle(i,:), 0];
            v1 = [quadrangle(mod(i-2,n)+1,:), 0];
            v2 = [quadrangle(mod(i,n)+1,:), 0];

            a = v1 - v2;
            b = pt - v2;
            distToLine = norm(cross(a,b)) / norm(a);

            if(distToLine <= minDist) 
                minDist = distToLine;
                removeIndex = i;
            end
        end
        quadrangle(removeIndex,:) = [];
    end
end

function [corners, middlePoints, middleImagePoints] = orderCorners(im, quadrangle)
    % ORDERCORNERS orders the 4 corner points of a quadrangle. The intention 
    % is to order the corners in a list as top left, top right, bottom right, bottom left.
    % In this implementation, the distance sum is computed between the ballot border 
    % middle points (points between corners) and the image border middle points (points
    % between image corners). The 4 corners points are orderd such that
    % this sum is minimal.
    % This seemed to be more stable than simply calculating the distances between ballot 
    % corners and image corners. Works well in 99.9% of cases.
    %
    % Author:
    %   Richard Binder
    %
    % Source:
    %   Self
    %
    % Inputs:
    %   im:     the image in which the simplifiedHull resides
    %   quadrangle:       a set of 4 corner points, making up a quadrangle
    %
    % Output:
    %   corners:                        the same set of points as in quadrangle, but ordered
    %   middlePoints:                the points between corners
    %   middleImagePoints:      the points between image corners

    middleImagePoints = [
    size(im, 2)/2, 1;
    size(im, 2), size(im, 1)/2;
    size(im, 2)/2, size(im, 1);
    1, size(im, 1)/2;
    ];     

    minDist = inf;
    corners = [];
    tmpCorners = quadrangle;
    middlePoints = quadrangle;
    for i = 1:4
        for j = 1:4
            tmpCorners(j,:) = quadrangle(mod(i+j-1,4)+1,:);
        end
        dist = 0;
        for j = 1:4
            middlePoints(j,:) = tmpCorners(j,:) + (tmpCorners(mod(j,4)+1,:) - tmpCorners(j,:))/2;
            dist = dist + norm(middlePoints(j,:) - middleImagePoints(j, :));
        end
        if(dist < minDist)
            minDist = dist;
            corners = tmpCorners;
            a=0;
        end
    end
end

function plotLinesAndCorners(im, lines, lines2, intersections, hull, corners, middlePoints, middleImagePoints)
    % PLOTLINESANDCORNERS plots all lines specified in lines and lines2,
    % and all points specified in intersections, hull, corners,
    % middlePoints, middleImagePoints on the image im
    %
    % Author:
    %   Richard Binder
    %
    % Source:
    %   Self
    
    global pltM;
    global pltN;
    global pltCount;

    subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
    hold on;
    imshow(im); title('Lines & Corners');
    plotLines(lines);
    plotLines(lines2);

    plot(intersections(:,1),intersections(:,2),'*', 'Color', 'blue');
    plot(hull(:,1),hull(:,2),'Color', 'red', 'lineWidth', 2.0);

    for j = 1:4
        plot(corners(j,1),corners(j,2),'o', 'Color', 'red', 'MarkerSize', 30);
        text(corners(j,1),corners(j,2), num2str(j),'Color', 'red','FontSize',20);
        
        plot(middlePoints(j,1),middlePoints(j,2),'o', 'Color', 'red', 'MarkerSize', 5);
        text(middlePoints(j,1),middlePoints(j,2), num2str(j),'Color', 'red','FontSize',5);
        plot(middleImagePoints(j,1),middleImagePoints(j,2),'o', 'Color', 'red', 'MarkerSize', 5);
        text(middleImagePoints(j,1),middleImagePoints(j,2), num2str(j),'Color', 'red','FontSize',5);
    end
    hold off;
end

function plotLines(lines)
    % PLOTLINES plots all lines specified in lines
    %
    % Author:
    %   Richard Binder
    %
    % Source:
    %   Self
    
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1), xy(:,2), 'LineWidth', 5, 'Color', 'black');
       plot(xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
       plot(xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'magenta');
    end
end


function rotated = correctTableRotation(ballotTable)
    % CORRECTTABLEROTATION corrects the rotation of ballotTable. 
    % It does this by first rotating it such that the horizontal length is
    % longer than the vertical length. Then it rotates it such that the
    % number of gradient Pixels is higher on the left side (since there is
    % more black pixels on the template of the ballot table there).
    %
    % Author:
    %   Richard Binder
    %
    % Source:
    %   Self
    %
    % Inputs: 
    %   ballotTable:    An image resembling a ballot table
    %
    % Output:
    %   rotated:        Same as ballotTable, but with corrected rotation
    
    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    global pltCount;
    
    rotated = ballotTable;
    
    grayTable = toGray(im2double(ballotTable));
    blurredTable = gaussfilt(grayTable, 5.0);
    % Create a gradient magnitude mask
    [gradMag, ~] = imgradient(blurredTable);
    gradThreshold = mean(gradMag(:)) - 0.2 * std(gradMag(:));
    maskedTable = (gradMag > gradThreshold);
    
    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
        imshow(maskedTable); title('Table Mask');
    end
    
    if size(maskedTable,1) > size(maskedTable,2) 
        maskedTable = rot90(maskedTable);
        rotated = rot90(rotated);
    end
    
    leftSide = maskedTable(:, 1:round(size(maskedTable,2)/2));
    rightSide = maskedTable(:, round(size(maskedTable,2)/2):size(maskedTable,2));
    if sum(leftSide(:)) < sum(rightSide(:))
        maskedTable = rot90(maskedTable);
        maskedTable = rot90(maskedTable);
        rotated = rot90(rotated);
        rotated = rot90(rotated);
    end
end