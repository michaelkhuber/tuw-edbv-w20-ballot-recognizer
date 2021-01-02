function transformed = Transform(im, ballotFilenameIn, step)
% TRANSFORM finds the ballot paper in the given input image and
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
%   im:        The input image
%
% Output:
%   transformed:    The transformed and cropped ballot

    global showPlot;
    global savePlot;
    global ballotFilename;
    global f;
    global pltM;
    global pltN;
    
    showPlot = true;
    savePlot = true;
    ballotFilename = ballotFilenameIn;
    pltM = 3;
    pltN = 5;
    pltCount = 1;

    if step ~= 1 && step ~= 2
        error("Invalid transform step");
    end

    %resize image to have 2000 pixels in height without deforming
    newSize = [2000,size(im,2)/size(im,1) * 2000];
    im = imresize(im,newSize);

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

    % normalize the image
    if step == 1
        [normalizedImage, pltCount] = Normalize(im, 20.0, true, pltCount);
    elseif step == 2
        [normalizedImage, pltCount] = Normalize(im, 5.0, false, pltCount);
    end

    % mask the image
    if step == 1
        [maskedImage, pltCount] = MaskImage(normalizedImage, pltCount);
    elseif step == 2
        [maskedImage, pltCount] = MaskImage2(normalizedImage, pltCount);
    end

    % get lines from hough transformation of the masked image
    [lines, lines2, pltCount] = HoughTransform(maskedImage, pltCount);

    % find intersections between vertical and horizontal lines
    intersections = findIntersections(im, lines, lines2);

    % find the convex hull from all the intersection points
    k = convhull(intersections);
    hull = [];
    hull(:,:) = intersections(k,:);

    % reduces the hull to 4 corner points
    simplifiedHull = findCorners(hull);

    % orders the corners as top-left, top-right, bottom-right, bottom-left
    [corners, middlePoints, middleImagePoints] = orderCorners(im, simplifiedHull);

    if(showPlot || savePlot) 
        pltCount = plotLinesAndCorners(im, lines, lines2, intersections, hull, corners, middlePoints, middleImagePoints, pltCount);
    end

    % Transforms the image such that the ballot corners make up a
    % (near) perfect rectangle.
    [imNew, pltCount] = HomographyTransformation(im, corners, pltCount);

    %resize ballot to have its original template size
    if step == 1
        newSize = [2500,3500];
    elseif step == 2
        newSize = [2150,3500];
    end

    imNew = imresize(imNew,newSize);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
        imshow(imNew); title('Resized & Cropped');
    end

    if(savePlot)
        print(f,sprintf("resources/results/Step2_Transform%i/Transform%i_%s", step, step, ballotFilename),'-dpng','-r700'); 
        clf(f);
        clear f;
    end

    transformed = imNew;
end





function intersections = findIntersections(im, lines, lines2)
    imH = size(im, 1);
    imW = size(im, 2);
    
    % Find the intersections of all the lines
    % Ignore the ones out-of-bound
    intersections = [];
    for i = 1:length(lines)
      for j = 1:length(lines2)
       %if( j <= i ), continue; end
        p1 = lines(i).rho;
        p2 = lines2(j).rho;
        t1 = lines(i).theta;
        t2 = lines2(j).theta;

        x = (p1*sind(t2)-p2*sind(t1))/(cosd(t1)*sind(t2)-sind(t1)*cosd(t2));
        y = (p1*cosd(t2)-p2*cosd(t1))/(sind(t1)*cosd(t2)-cosd(t1)*sind(t2));
        if (isnan(x) || isnan(y) || x < 1 || x > imW || y < 1 || y > imH)
            continue;
        end
        intersections = [intersections; x, y];
      end
    end
end

function simplifiedHull = findCorners(hull)
    %Reduce the points inside the hull set until there is only 4 points
    %this should be the corner points of our rectangle hull
    simplifiedHull = hull;
    while (length(simplifiedHull) > 4)
        minDist = inf;
        n = length(simplifiedHull);
        removeIndex = 0;

        for i=1:n
            pt = [simplifiedHull(i,:), 0];
            v1 = [simplifiedHull(mod(i-2,n)+1,:), 0];
            v2 = [simplifiedHull(mod(i,n)+1,:), 0];

            a = v1 - v2;
            b = pt - v2;
            distToLine = norm(cross(a,b)) / norm(a);

            if(distToLine <= minDist) 
                minDist = distToLine;
                removeIndex = i;
            end
        end
        simplifiedHull(removeIndex,:) = [];
    end
end

function [corners, middlePoints, middleImagePoints] = orderCorners(im, simplifiedHull)
    % Order the 4 ballot corner points such that the sum of the distance to the outter
    % image borders is as small as possible. The intention is to order the corners 
    % in a list as top left, top right, bottom right, bottom left.
    % In this implementation, the distance sum is computed between the ballot border middle points 
    % (points between corners) and the image border middle points. This seemed to be more stable 
    % than simply calculating the distances between ballot corners and image corners. Works well in 99.9%
    % of cases.

    middleImagePoints = [
    size(im, 2)/2, 1;
    size(im, 2), size(im, 1)/2;
    size(im, 2)/2, size(im, 1);
    1, size(im, 1)/2;
    ];     

    minDist = inf;
    corners = [];
    tmpCorners = simplifiedHull;
    middlePoints = simplifiedHull;
    for i = 1:4
        for j = 1:4
            tmpCorners(j,:) = simplifiedHull(mod(i+j-1,4)+1,:);
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

function pltCount = plotLinesAndCorners(im, lines, lines2, intersections, hull, corners, middlePoints, middleImagePoints, pltCount)
    global pltM;
    global pltN;

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
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1), xy(:,2), 'LineWidth', 5, 'Color', 'black');
       plot(xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
       plot(xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'magenta');
    end
end