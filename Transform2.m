% AUTHORS
% Jakob @ TUWIEN
% Binder Richard @ TUWIEN
function transformed = Transform2(im, ballotFilenameIn)
        global showPlot;
        global savePlot;
        global ballotFilename;
        global f;
        global pltM;
        global pltN;
        
        showPlot = false;
        savePlot = false;
        ballotFilename = ballotFilenameIn;
        pltM = 3;
        pltN = 5;
        pltCount = 1;
        
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
        
        [normalizedImage, pltCount] = normalize(im, pltCount);
        
        minPeaks = 10;
        maxPeaks = 50;
        
        [maskedImage, pltCount] = maskImage(normalizedImage, pltCount);
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
        
        
        %[intersections, pltCount] = removeIsolated(intersections, maskedImage, pltCount);
        

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
        
        %resize ballot to have its original template size dimensions of 2500x3500
        %alternative: keep dimensions of input image
        %newSize = [2000,size(imNew,2)/size(imNew,1) * 2000];
        newSize = [2150,3500];
        imNew = imresize(imNew,newSize);

        if(showPlot || savePlot) 
            % Plot the results
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(imNew); title('Resized & Cropped');
        end
        
        if(savePlot)
            print(f,strcat("resources/results/Transform2_", ballotFilename),'-dpng','-r700'); 
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
        componentMask = gradMask;
		cc = bwconncomp(gradMask, 4);
        
        ccSizes = zeros(cc.NumObjects,1);
        for i = 1 : cc.NumObjects
            currCC = cc.PixelIdxList{i};
            ccSizes(i) = size(currCC, 1);
        end
        
		ccSizeThreshold = mean(ccSizes(:)) + strength*std(ccSizes(:));
        maxComponent = [];
        
        for i = 1 : cc.NumObjects
            currCC = cc.PixelIdxList{i};
            if size(currCC, 1) > size(maxComponent, 1)
                maxComponent = currCC;
            end
            if size(currCC, 1) < ccSizeThreshold
                componentMask(currCC) = 0;
            end
        end
        
        componentMask(maxComponent) = 1;
end

function [maskedImage, pltCount] = maskImage(img, pltCount)
        global showPlot;
        global savePlot;
        global pltM;
        global pltN;
       

		% Create a gradient magnitude mask
		[gradMag, ~] = imgradient(img);
		gradThreshold = mean(gradMag(:)) - 0.2 * std(gradMag(:));
		gradMask = (gradMag > gradThreshold);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(gradMask); title('Grad Mask');
        end
        
		% Find all the connected compoments & remove small ones
        componentMask = reduceComponents(gradMask, 3.0);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
            imshow(componentMask); title('Component Reduction');
        end

        DilationMask = componentMask;
        se = strel('octagon',3);
        DilationMask = imdilate(DilationMask, se);
%         se = strel('octagon',3);
%         DilationMask = imerode(DilationMask, se);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(DilationMask); title('Dilation Mask');
        end        
        
        ConvexHull = bwconvhull(DilationMask);
        ConvexHull(:,size(ConvexHull,2)) = 0.0;
        ConvexHull(:,1) = 0.0;
        ConvexHull(size(ConvexHull,1),:) = 0.0;
        ConvexHull(1, :) = 0.0;
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(ConvexHull); title('Convex Hull');
        end
        
        %ConvexHull = imgaussfilt(ConvexHull, 20);
		[gradMag, ~] = imgradient(ConvexHull);
        ConvexHull = gradMag > 0.00001;
        
        se = strel('octagon',18);
        ConvexHull = imdilate(ConvexHull, se);
        
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(ConvexHull); title('Border Convex Hull');
        end
        
        
        
        maskedImage = ConvexHull;
        
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

function [normalizedImage, pltCount] = normalize(image, pltCount)
    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    
    [nobg, saturationSuccess, pltCount] = removeBackground(image, pltCount, 0.15, 0.01, false, "Saturation Histogram Analysis");
    if ( saturationSuccess )        
        nobg = im2double(nobg);
        nobg = rgb2gray(nobg);
        blurredImg = imgaussfilt(nobg, 5.0);
        normalizedImage = blurredImg;
        return;
    end
        
    [nobg, saturationSuccess, pltCount] = removeBackground(image, pltCount, 0.15, 0.3, true, "Extended Saturation Histogram Analysis");
    if ( saturationSuccess )
        nobg = im2double(nobg);
        nobg = rgb2gray(nobg);    
        
        chroma = getChroma(image);
        white = sqrt(3*double(255)^2);
        chroma = sqrt(double(chroma(:,:,1)).^2 + double(chroma(:,:,2)).^2 + double(chroma(:,:,3)).^2) ./ white;
        chroma = adapthisteq(chroma,'clipLimit',0.02,'Distribution','rayleigh');
        
        chroma(nobg < 0.1) = chroma(nobg < 0.1) - 3 * std(chroma(:));
        chroma( chroma < 0.01 ) = 0.0;
        gray = rgb2gray(im2double(image));
        %gray = sqrt( chroma .* gray );
        gray(nobg < 0.1) = gray(nobg < 0.1) - 2 * std(gray(:));
        blurredImg = imgaussfilt(gray, 100);
        blurredImg(nobg > 0.1) = gray(nobg > 0.1);
        blurredImg = imgaussfilt(blurredImg, 5.0);

        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(blurredImg); title("Background Based Chromaticity Blurring");
        end
        normalizedImage = blurredImg;
        return;
    end
    
    grayImg = rgb2gray(im2double(image));
    blurredImg = imgaussfilt(grayImg, 5.0);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        imshow(blurredImg); title("Fallback To Denoised Brightness");
    end
    
    normalizedImage = blurredImg;
end

%Removes Background based on the histogram of the saturation.
%On a local minima of the saturation histogram, a threshold is applied onto
%the image.
function [thresholdImage, success, pltCount] = removeBackground(image, pltCount, maxProminence, minProminence, doZeroApprox, histogramTitle)
    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    
    hsv = rgb2hsv(image);    
    mask = 1.0 - hsv(:,:,2);
    
    [counts, edges] = histcounts(mask, 100);
    [locMin, locMax] = getLocals(counts, maxProminence, minProminence);
    edges = edges(1:(length(edges)-1));
    
    if( isempty(locMin) ) 
        thresholdImage = image;
        success = 0;
        return;
    end    
    
    localMinDist = edges(locMin);
    
    if( doZeroApprox )
        zeroApprox = locMin - 5 * counts(locMin) * (locMax - locMin) / (counts(locMax) - counts(locMin));
        zeroApprox = round( zeroApprox );
        if( zeroApprox < 1 )
            zeroApprox = 1;
        end
        localMinDist = edges(round(zeroApprox));
    end
    
    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        hold on;
        plot(1.0-edges,counts,'Color', 'red');
        plot(1.0-edges(locMax),counts(locMax),'^','Color', 'blue');
        plot(1.0-edges(locMin),counts(locMin),'v','Color', 'blue');
        
        title(histogramTitle);
        hold off;
    end
    
    binary = mask <= localMinDist;
    pixelsRemoved = sum(binary(:));
    
    thresholdImage = image;    
    for i = 1:3
        channel = image(:,:,i);
        channel(binary) = 0.0;
        thresholdImage(:,:,i) = channel;
    end
    
    % assume that number of removed pixels is an indicator whether the
    % background removal succeeded or not
    minRemoved = 0.01 * length(binary(:));
    maxRemoved = 0.15 * length(binary(:));
    
    if((pixelsRemoved < minRemoved) || (pixelsRemoved > maxRemoved) )
        success = 0;
    else
        success = 1;
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(thresholdImage); title("Saturation Based Background Removal");
        end
    end
end

%computes a local minimum and local maximum with difference greater than
%prominence
function [locMin, locMax] = getLocals(counts, maxProminence, minProminence)
    maxProminence = maxProminence * max(counts(:));

    %get maxima
    maxima = [];
    for i = 2:(length(counts)-1)
        if( (counts(i) > counts(i-1)) && (counts(i) > counts(i+1)) )
            maxima = [maxima, i];
        end
    end
    
    %edge cases
    if( counts(1) > counts(2) )
        maxima = [1, maxima];
    end
    if( counts(length(counts)) > counts(length(counts)-1) )
        maxima = [maxima, length(counts)];
    end
    
    minima = [];
    for i = 2:(length(counts)-1)
        if( (counts(i) < counts(i-1)) && (counts(i) < counts(i+1)) )
            minima = [minima, i];
        end
    end
    
    %edge cases
    if( counts(1) < counts(2) )
        minima = [1, minima];
    end
    if( counts(length(counts)) < counts(length(counts)-1) )
        minima = [minima, length(counts)];
    end
    
    prominentMaxima = [];
    prominentMinima = [];
    for i = 1:length(maxima)
        if( counts(maxima(i)) > maxProminence )
            for j = 1:length(minima)
                if( minima(j) < maxima(i) && counts(minima(j)) < minProminence * counts(maxima(i)) )
                    prominentMaxima = [prominentMaxima, maxima(i)];
                    prominentMinima = [prominentMinima, minima(j)];
                end
            end
        end
    end
    
    %get right most prominent maximum
    locMax = max(prominentMaxima);
    
    %and its right most assigned prominent minimum
    locMin = max(prominentMinima(prominentMaxima==locMax));
end

%gets Chromaticity of a (double) rgb image
function chroma = getChroma(image)
    image = im2double(image);
    brightness = rgb2gray(image);
    brightness(brightness < 0.01) = 0.01;
    r = image(:,:,1) ./ brightness;
    g = image(:,:,2) ./ brightness;
    b = image(:,:,3) ./ brightness;
    chroma = cat(3, r, g, b);
    chroma = im2uint8(chroma);
end