% AUTHORS
% Jakob @ TUWIEN
% Binder Richard @ TUWIEN
function transformed = Transform(im, ballotFilenameIn)
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
        pltN = 4;
        pltCount = 1;
        
        %resize image to have 2000 pixels in height without deforming
        newSize = [2000,size(im,2)/size(im,1) * 2000];
        im = imresize(im,newSize);
            
        if(showPlot)
            f = figure(1);
            clf('reset');
        elseif(savePlot)
            f = figure('visible','off');
        end
        
        if(showPlot || savePlot)
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(im); title('Original');
        end
        
        [noBgrdImg, pltCount] = cDarkendBackground(im, pltCount);
		[imH, imW] = size(noBgrdImg);
        
%         % find an image mask in which there is a given number of
%         % components
%         % this includes tuning parameters for applying a gauss filter,
%         % applying a median filter, peforming component reduction and
%         % performing dilation onto vertical/horizontal lines
%         % these parameters are tuned until we find given component numbers
%         % the implementation is improvable, so it might not always work 100% well
%         minComponents = 500;
%         maxComponents = 10000;
%         med = -2;
%         numComponents = inf;
%         while(numComponents < minComponents || numComponents > maxComponents)
%             med = med+2;
%             [maskedImage, numComponents, pltCount2] = maskImage(noBgrdImg, pltCount, 5, med, 0, 1);
%             if(med>=8) 
%                 break;
%             end
%         end
%         
%         minComponents = 50;
%         maxComponents = 100;
%         dil = -1;
%         numComponents = inf;
%         while(numComponents > maxComponents)
%             dil = dil+2;
%             [maskedImage, numComponents, pltCount2] = maskImage(noBgrdImg, pltCount, 5, med, 0, dil);
%             if(dil>=6) 
%                 break;
%             end
%         end
%         
%         maxComponents = 60;
%         sd = -2;
%         numComponents = inf;
%         while(numComponents > maxComponents)
%             sd = sd+1;
%             [maskedImage, numComponents, pltCount2] = maskImage(noBgrdImg, pltCount, 5, med, sd, dil);
%             if(sd>=3) 
%                 break;
%             end
%         end
%         
%         % perform hough transformation onto the masked image
%         % if we cannot find enough peaks, then we again tune image mask
%         % parameters until we find enough peaks
%         % this is also improvable, and might not always work 100% well
%         minPeaks = 10;
%         maxPeaks = 30;
%         [H, theta, rho, H2, theta2, rho2]  = houghTransform(maskedImage, 15);
%         P = houghpeaks(H, maxPeaks);
%         P2 = houghpeaks(H2, maxPeaks);
%         if(length(P) < minPeaks || length(P2) < minPeaks)
%             [maskedImage, numComponents, pltCount2] = maskImage(noBgrdImg, pltCount, 5, med, sd, dil+2);
%             [H, theta, rho, H2, theta2, rho2]  = houghTransform(maskedImage, 15);
%             P = houghpeaks(H, maxPeaks);
%             P2 = houghpeaks(H2, maxPeaks);
%         end
        
        minPeaks = 10;
        maxPeaks = 20;
        
        
        [maskedImage, numComponents, pltCount2] = maskImage2(noBgrdImg, pltCount, 5, 0, 0, 0);
        [H, theta, rho, H2, theta2, rho2]  = houghTransform(maskedImage, 5);
        P = houghpeaks(H, maxPeaks);
        P2 = houghpeaks(H2, maxPeaks);
        
        
        
        %set pltCount to new value
        pltCount = pltCount2;
        
        
        
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
		lines = houghlines(maskedImage, theta, rho,P, 'FillGap', 200, 'MinLength', 1000);
		lines2 = houghlines(maskedImage, theta2, rho2,P2, 'FillGap', 200, 'MinLength', 1000);
        
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
			if (isnan(x) || isnan(y) || x <= 0 || x > imW || y <= 0 || y > imH)
                continue;
            end
			intersections = [intersections; x, y];
		  end
        end
        
        
        [intersections, pltCount] = removeIsolated(intersections, pltCount);
        

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
        
        % Order the 4 corner points so that the sum of the distance to the outter
        % image borders is as small as possible
        % the intention is to order the corners as top left, top right,
        % bottom right, bottom left
        % currently, this just calculates the points in the middle of the
        % lines between the corners, and then calculates the sum of the
        % distances to the image border middle points
        % this seemed to be more stable than simply calculating the
        % distances between corners and image corners
        % might be improvable, but works well in ~95% of cases
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

		% Measure the skewed widths & heights
		heightL = norm(corners(1,:) - corners(4,:));
		heightR = norm(corners(2,:) - corners(3,:));
		widthT = norm(corners(1,:) - corners(2,:));
		widthB = norm(corners(3,:) - corners(4,:));

		% Set up the target image dimensions
		% Use the maximum of skewed width and height 
		% to approxmate the target dimensions
		imNewHeight = max([heightL, heightR]);
		imNewWidth  = max([widthT, widthB]);
        
		cornersNew = [         1,           1; 
					  imNewWidth,           1;
					  imNewWidth, imNewHeight;
							   1, imNewHeight];

		% Compute the homography matrix
		corners = corners';
		cornersNew = cornersNew';
                
		h = ComputeHNorm(cornersNew, corners);
        
		% Apply it to the original image
		tform = projective2d(h');
		[imNew, RB] = imwarp(im, tform);
        
        if(showPlot || savePlot) 
            % Plot the results
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(imNew); title('Transformed');
        end
        
        %Crop the image to to the 4 corners of the ballot
        upper_left  = corners([1,2],1)';
        lower_right = corners([1,2],3)';
        
        [lower_right_X,lower_right_Y] = transformPointsForward(tform, lower_right(1), lower_right(2));
        [X1, Y1] = worldToIntrinsic(RB, lower_right_X, lower_right_Y);
        
        [upper_left_X, upper_left_Y] = transformPointsForward(tform, upper_left(1), upper_left(2));
        [X2, Y2] = worldToIntrinsic(RB, upper_left_X, upper_left_Y);
        
        imNew = imcrop(imNew, [X2 Y2 X1-X2 Y1-Y2]);
        
        %resize image to have 2000 pixels in height without deforming
        newSize = [2000,size(imNew,2)/size(imNew,1) * 2000];
        imNew = imresize(imNew,newSize);

        if(showPlot || savePlot) 
            % Plot the results
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(imNew); title('Cropped');
        end
        
        if(savePlot)
            print(f,strcat("resources/results/Transform_", ballotFilename),'-dpng','-r700'); 
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

function H2to1 = ComputeHNorm(p1, p2)
        % compute the homography norm
		[row,col] = size(p2); 
		L = [0]; z13 = [0,0,0];
		p2 = [p2;ones(1,col)];
		p1 = [p1;ones(1,col)];
		t1 = norm_matrix(p2);
		t2 = norm_matrix(p1);
		p2 = t1*p2; 
		p1 = t2*p1;

		L = [p2(:,1).', z13, p1(1,1)*(-1)*(p2(:,1).');
			 z13, p2(:,1).', p1(2,1)*(-1)*(p2(:,1).')];
		for i=2:col
			L = [L;
				 p2(:,i).', z13, p1(1,i)*(-1)*(p2(:,i).');
				 z13, p2(:,i).', p1(2,i)*(-1)*(p2(:,i).')];
		end

		a = (L.')*L;
		[v,d] = eig(a);

		min = d(9,9);
		index = 0;
		for i=1:9
			if (d(i,i)~=0) & (d(i,i)<min)
				min = d(i,i);
				index = i;
			end
		end

		x = v(:,index).';

		Hnorm = [x(1),x(2),x(3);x(4),x(5),x(6);x(7),x(8),x(9)];
		H2to1 = inv(t2)*Hnorm*t1;%turn "normalized H" back to H
		end

function Tnorm = norm_matrix(p2)
		[row,col] = size(p2); %row=2, col=n
		avg = sum(p2,2)/col;
		totaldist = 0;
		for i=1:col
			totaldist = totaldist + pdist([avg(1),avg(2); p2(1,i),p2(2,i)]);
		end
		k = 1/(totaldist/(col*sqrt(2))); 
		Ttran =  [1 0 -avg(1);0 1 -avg(2);0 0 1];
		Tscale = [k 0 0;0 k 0;0 0 1];
		Tnorm = Tscale*Ttran;
end

%reduce components of input mask
% numComponents is the number of remaining components
function componentMask = reduceComponents(gradMask, strength)
		cc = bwconncomp(gradMask, 4);
        
        ccSizes = zeros(cc.NumObjects,1);
        for i = 1 : cc.NumObjects
		  currCC = cc.PixelIdxList{i};
		  ccSizes(i) = size(currCC, 1);
        end
        
		ccSizeThreshold = mean(ccSizes(:)) + strength*std(ccSizes(:));
        
		for i = 1 : cc.NumObjects
		  currCC = cc.PixelIdxList{i};
		  if size(currCC, 1) < ccSizeThreshold
			gradMask(currCC) = 0;
		  end
        end
        componentMask = gradMask;
end

function [newImage, pltCount] = cDarkendBackground(image, pltCount)
    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    
    grayImg = rgb2gray(image);
    grayImg = im2double(grayImg);
    
    % equalize/normalize the image to get a better distribution
    normalized = adapthisteq(grayImg,'clipLimit',0.02,'Distribution','rayleigh');
    
    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        imshow(normalized); title("Histogram Equalization");
    end
    
    newImage = normalized;
end



function [maskedImage, numComponents, pltCount] = maskImage2(img, pltCount, gaussStrength, medianStrength, componentReducingStrength, dilationStrength)
        global showPlot;
        global savePlot;
        global pltM;
        global pltN;
        
		% Blur the image for denoising
		H = fspecial('gaussian', [20, 20], 20);
		blurredImg = imfilter(img, H, 'replicate');
        if(medianStrength > 0)
            blurredImg = medfilt2(blurredImg, [medianStrength, medianStrength]);
        end
        blurredImg = medfilt2(blurredImg, [30, 30]);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(blurredImg); title('Denoising Filters');
        end

		% Create a gradient magnitude mask
		[gradMag, ~] = imgradient(blurredImg);
		gradThreshold = mean(gradMag(:)) + std(gradMag(:));
		gradMask = (gradMag > gradThreshold);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(gradMask); title('Grad Mask');
        end
        
		% Find all the connected compoments & remove small ones
        componentMask = reduceComponents(gradMask, 2);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(componentMask); title('Component Reduction');
        end

        
% 		% Find edge mask
%         edgeMask = edge(gradMask, 'canny');
        
		DilationMask = imfilter(componentMask, H, 'replicate');
        se = strel('line',100,100);
        se2 = strel('line',100,100);
        DilationMask = imdilate(DilationMask, se);
        %DilationMask = imerode(DilationMask, se2);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(DilationMask); title('Dilation Mask');
        end
        
        maskedImage = DilationMask;
		cc = bwconncomp(maskedImage, 4);
        numComponents = cc.NumObjects;
end

function [maskedImage, numComponents, pltCount] = maskImage(img, pltCount, gaussStrength, medianStrength, componentReducingStrength, dilationStrength)
        global showPlot;
        global savePlot;
        global pltM;
        global pltN;
        
		% Blur the image for denoising
		H = fspecial('gaussian', [gaussStrength, gaussStrength], gaussStrength);
		blurredImg = imfilter(img, H, 'replicate');
        if(medianStrength > 0)
            blurredImg = medfilt2(blurredImg, [medianStrength, medianStrength]);
        end
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(blurredImg); title('Denoising Filters');
        end

		% Create a gradient magnitude mask
		gradThreshold = mean(blurredImg(:)) - std(blurredImg(:));
		[gradMag, ~] = imgradient(blurredImg);
		gradMask = (gradMag > gradThreshold);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(gradMask); title('Grad Mask');
        end
        
		% Find all the connected compoments & remove small ones
        componentMask = reduceComponents(gradMask, componentReducingStrength);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(componentMask); title('Component Reduction');
        end

        
% 		% Find edge mask
%         edgeMask = edge(gradMask, 'canny');
        
		DilationMask = imfilter(componentMask, H, 'replicate');
        se = strel('line',dilationStrength,dilationStrength*4);
        DilationMask = imdilate(componentMask, se);
        se = strel('line',dilationStrength*4,dilationStrength);
        DilationMask = imdilate(DilationMask, se);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(DilationMask); title('DilationEdge Mask');
        end
        
        maskedImage = DilationMask;
		cc = bwconncomp(maskedImage, 4);
        numComponents = cc.NumObjects;
end

function [H, theta, rho, H2, theta2, rho2] = houghTransform(maskedImg, precision)
		% Use Hough transform to find vertical borders: Theta goes from 0°
		% to 90°
		[H, theta, rho] = hough(maskedImg, 'RhoResolution', precision, 'Theta', [0.0:0.5:85.5]);
        
		% Use Hough transform to find horizontal borders: Theta goes from -90°
		% to 0°
		[H2, theta2, rho2] = hough(maskedImg, 'RhoResolution', precision, 'Theta', [-90.0:0.5:0.0]);
end


function [newPoints, pltCount] = removeIsolated(points, pltCount)
    global showPlot;
    global savePlot;
    global pltN;
    global pltM;

    minX = min(points(:,1));
    minY = min(points(:,2));
    maxX = max(points(:,1));
    maxY = max(points(:,2));
    
    % mid = [mean(points(:,1)), mean(points(:,2))];
    mid = [(maxX-minX)/2, (maxY-minY)/2];
    corners = [minX maxY; maxX maxY; maxX minY; minX minY];
    
    if(showPlot || savePlot)
        subplot(pltM, pltN, pltCount);
        cla();
        hold on;
        plot(points(:,1), points(:,2),'*', 'Color', 'blue'); title("Removing isolated Points");
        plot(mid(1), mid(2),'*', 'Color','red');
        hold off;
    end
    
    distance = vecnorm((points - mid)')';
    sd = std(distance);
    m = mean(distance);
    
    newPoints = points(distance < m+1*sd,:);
    
    best = newPoints;
    last = newPoints;
    
    k = convhull(last);
    hull = last(k,:);
    lastArea = polyarea(hull(:,1),hull(:,2));
    bestArea = inf;
    initialArea = lastArea;
    
    relation = 0;
    while(relation < 0.95 && lastArea/initialArea > 0.7)
        if(showPlot || savePlot)
            subplot(pltM, pltN, pltCount);
            hold on;
            plot(hull(:,1), hull(:,2),'Color', 'red');
            hold off;
        end
        best = last;
        bestArea = lastArea;
        
        last(k,:) = [];
        k = convhull(last);
        hull = last(k,:);
        lastArea = polyarea(hull(:,1),hull(:,2));
        relation = lastArea/bestArea;
    end
    
    if(showPlot || savePlot)
         pltCount = pltCount + 1;
    end
    
    newPoints = best;
end

