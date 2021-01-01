% AUTHORS
% Jakob @ TUWIEN
% Binder Richard @ TUWIEN
function transformed = Transform(im, ballotFilenameIn, step)
% TRANSFORM determines the ballot paper in the given input image,
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
        
        showPlot = false;
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
		imH = size(im, 1);
		imW = size(im, 2);
            
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
        
        %perform hough transformation
        maxPeaks = 50;
        
        [H, theta, rho, H2, theta2, rho2]  = houghTransform(maskedImage, 5);
        P = houghpeaks(H, maxPeaks);
        P2 = houghpeaks(H2, maxPeaks);       
        
        
        % combine horizontal and vertical lines
        HCombined = [H2,H];
        thetaCombined = [theta2, theta];
        rhoCombined = [rho2, rho];
        PCombined = [P2;[P(:,1),P(:,2)+size(H2,2)]];
        
        if(showPlot || savePlot)
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(imadjust(rescale(HCombined)),'XData',thetaCombined,'YData',rhoCombined,...
          'InitialMagnification','fit'); %plot Hough transformation
             title('Hough Transform');
             xlabel('\theta'), ylabel('\rho');
             axis on, axis normal, hold on;
            plot(thetaCombined(PCombined(:,2)),rhoCombined(PCombined(:,1)),'s','color','white'); %plot peaks
            colormap(gca,hot);
        end
        
        % Get Lines from Hough Transformation
		lines = houghlines(maskedImage, theta, rho,P, 'FillGap', 500, 'MinLength', 500);
		lines2 = houghlines(maskedImage, theta2, rho2,P2, 'FillGap', 500, 'MinLength', 500);
        
        lines = [lines, lines2];

		% Find the intersections of all the lines
		% Ignore the ones out-of-bound
		intersections = [];
		for i = 1:length(lines)
		  for j = 1:length(lines)
            if( j <= i ), continue; end
			p1 = lines(i).rho;
			p2 = lines(j).rho;
			t1 = lines(i).theta;
			t2 = lines(j).theta;

			x = (p1*sind(t2)-p2*sind(t1))/(cosd(t1)*sind(t2)-sind(t1)*cosd(t2));
			y = (p1*cosd(t2)-p2*cosd(t1))/(sind(t1)*cosd(t2)-cosd(t1)*sind(t2));
			if (isnan(x) || isnan(y) || x < 1 || x > imW || y < 1 || y > imH)
                continue;
            end
			intersections = [intersections; x, y];
		  end
        end
        
        if step == 1
            [intersections, pltCount] = removeIsolated(intersections, maskedImage, pltCount);
        end
        

        % find the convex hull from all the intersection points
		k = convhull(intersections);
        hull = [];
        hull(:,:) = intersections(k,:);
        
        
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
        
        middleImagePoints = [
            size(im, 2)/2, 1;
            size(im, 2), size(im, 1)/2;
            size(im, 2)/2, size(im, 1);
            1, size(im, 1)/2;
            ];     
        
        % Order the 4 ballot corner points such that the sum of the distance to the outter
        % image borders is as small as possible. The intention is to order the corners 
        % in a list as top left, top right, bottom right, bottom left.
        % In this implementation, the distance sum is computed between the ballot border middle points 
        % (points between corners) and the image border middle points. This seemed to be more stable 
        % than simply calculating the distances between ballot corners and image corners. Works well in 99.9%
        % of cases.
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
        
                
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            hold on;
            imshow(im); title('Lines & Corners');
            plotLines(lines);
            plotLines(lines2);
            
            plot(intersections(:,1),intersections(:,2),'*', 'Color', 'blue');
            plot(hull(:,1),hull(:,2),'Color', 'red', 'lineWidth', 2.0);
            
            for j = 1:4
                plot(corners(j,1),corners(j,2),'o', 'Color', 'red', 'MarkerSize', 30);
                text(corners(j,1),corners(j,2), num2str(j),'Color', 'blue','FontSize',20);
%                 middlePoints(j,:) = corners(j,:) + (corners(mod(j,4)+1,:) - corners(j,:))/2;
%                 plot(middleCornerPoints(j,1),middleCornerPoints(j,2),'o', 'Color', 'red', 'MarkerSize', 30);
%                 text(middleCornerPoints(j,1),middleCornerPoints(j,2), num2str(j),'Color', 'red','FontSize',30);
%                 plot(middleImagePoints(j,1),middleImagePoints(j,2),'o', 'Color', 'red', 'MarkerSize', 30);
%                 text(middleImagePoints(j,1),middleImagePoints(j,2), num2str(j),'Color', 'red','FontSize',30);
            end
            hold off;
        end
       
        % Transforms the image such that the ballot corners make up a
        % (near) perfect rectangle.
        [imNew, pltCount] = HomographyTransformation(im, corners, pltCount);
        
        %resize ballot to have its original template size dimensions of 2500x3500
        if step == 1
            newSize = [2500,3500];
        elseif step == 2
            newSize = [2150,3500];
        end
        
        imNew = imresize(imNew,newSize);

        if(showPlot || savePlot) 
            % Plot the results
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

function plotLines(lines)
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'green');
       plot(xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
       %plot(xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'red');
    end
end

function [H, theta, rho, H2, theta2, rho2] = houghTransform(maskedImg, precision)
		% Use Hough transform to find vertical borders: Theta goes from 0°
		% to 90°
		[H, theta, rho] = hough(maskedImg, 'RhoResolution', precision, 'Theta', [0.0:0.3:85.5]);
        
		% Use Hough transform to find horizontal borders: Theta goes from -90°
		% to 0°
		[H2, theta2, rho2] = hough(maskedImg, 'RhoResolution', precision, 'Theta', [-90.0:0.3:0.0]);
end

function [newPoints, pltCount] = removeIsolated(points, mask, pltCount)
    global showPlot;
    global savePlot;
    global pltN;
    global pltM;
    
    se = strel('line',50,0);
    DilationMask1 = imdilate(mask, se);
    se = strel('line',50,90);
    DilationMask2 = imdilate(mask, se);
    mask = DilationMask1 | DilationMask2;

    minX = min(points(:,1));
    minY = min(points(:,2));
    maxX = max(points(:,1));
    maxY = max(points(:,2));
    
    % mid = [mean(points(:,1)), mean(points(:,2))];
    % mid = [(maxX-minX)/2, (maxY-minY)/2];
    
    if(showPlot || savePlot)
        subplot(pltM, pltN, pltCount);
        hold on;
        cla();
        plot(points(:,1), points(:,2),'.', 'Color', 'blue'); title("Removing isolated Points");
        set(gca, 'YDir','reverse');
        hold off;
    end
    
    % make sure that the points dont clip outside the borders
    maskedPoints = round(points);
    maskedPoints(maskedPoints <= 1) = 1;
    maskedPoints(maskedPoints(:,1)>=size(mask,2)-1,1) = size(mask,2);
    maskedPoints(maskedPoints(:,2)>=size(mask,1)-1,2) = size(mask,1);
    
    % define black image
    ptsMask = mask;
    ptsMask(:,:) = 0;
    % get linear indices for all the points
    linearIndices = sub2ind(size(mask), maskedPoints(:,2), maskedPoints(:,1));
    % set all pixels at the linearIndices the same as the mask
    ptsMask(linearIndices) = mask(linearIndices);
    %get rows and columns of all pixels that are white ( greater than 0.9
    %is just cause double comparison doesnt work )
    [r, c] = find(ptsMask > 0.9);
    %now we have all points that originally were at a white pixel on the
    %mask
    maskedPoints = [c, r];
    newPoints = maskedPoints;
    
    if(showPlot || savePlot)
        subplot(pltM, pltN, pltCount);
        hold on;
        plot(maskedPoints(:,1), maskedPoints(:,2),'*', 'Color', 'green');
        hold off;
    end
    
%     mid = [mean(newPoints(:,1)), mean(newPoints(:,2))];
%     distance = vecnorm((newPoints - mid)')';
%     sd = std(distance);
%     m = mean(distance);
%     closePoints = newPoints(distance < m+3*sd,:);
%     newPoints = closePoints;
%     
    k = convhull(newPoints);
    hull = newPoints(k,:);
    
    if(showPlot || savePlot)
        subplot(pltM, pltN, pltCount);
        hold on;
        plot(hull(:,1), hull(:,2),'LineWidth',1,'Color', 'red','LineStyle','-');
        hold off;
    end
    
    if(showPlot || savePlot)
         pltCount = pltCount + 1;
    end
    
    newPoints = newPoints;
end
