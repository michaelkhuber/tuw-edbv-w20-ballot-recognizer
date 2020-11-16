function transformed = Transform(im)
        imGray = rgb2gray(im);
		imGray = im2double(imGray);
		[imH, imW] = size(imGray);

		% Blur the image for denoising
		H = fspecial('gaussian', [5, 5], 5);
		imGray = imfilter(imGray, H, 'replicate');

		% Create a gradient magnitude mask
		gradThreshold = 0.3;
		[gradMag, ~] = imgradient(imGray);
		gradMask = (gradMag > gradThreshold);

		% Find all the connected compoments & remove small ones
		cc = bwconncomp(gradMask);
		ccSizeThreshold = 1000;
		for i = 1 : cc.NumObjects
		  currCC = cc.PixelIdxList{i};
		  if size(currCC, 1) < ccSizeThreshold
			gradMask(currCC) = 0;
		  end
		end

		% Find the mask for foreground with convex hull
		foregroundMask = bwconvhull(gradMask);
		edgeMask = edge(foregroundMask, 'Sobel');

		% Use Hough transform to find the borders of the card
		[H, theta, rho] = hough(edgeMask, 'RhoResolution', 5, 'Theta', [-90:0.5:89.5]);
		P = houghpeaks(H, 100); % 100 is just an arbitrarily large number
		lines = houghlines(edgeMask, theta, rho,P, 'FillGap', 30, 'MinLength', 3);

		% Find the intersections of all the lines
		% Ignore the ones out-of-bound
		corners = [];
		for i = 1:length(lines)
		  for j = 1:length(lines)
			if i>=j, continue; end;
			p1 = lines(i).rho;
			p2 = lines(j).rho;
			t1 = lines(i).theta;
			t2 = lines(j).theta;

			x = (p1*sind(t2)-p2*sind(t1))/(cosd(t1)*sind(t2)-sind(t1)*cosd(t2));
			y = (p1*cosd(t2)-p2*cosd(t1))/(sind(t1)*cosd(t2)-cosd(t1)*sind(t2));
			if x <= 0 || x > imW || y <= 0 || y > imH, continue; end;
			corners = [corners; x, y];
		  end
		end

		% Re-order corners this way: tl, tr, br, bl
		% Assume that the tl corner is closest to 1,1, etc.
		imageCorners = [          1,           1;
						size(im, 2),           1;
						size(im, 2), size(im, 1);
								  1, size(im, 1)];
		cornersTmp = [];
		for i = 1 : 4
		  cornersVector = corners - repmat(imageCorners(i, :), size(corners, 1), 1);
		  dist = (cornersVector(:, 1).^2 + cornersVector(:, 2).^2) .^ 0.5;
		  [~, ind] = min(dist);
		  cornersTmp(i, :) = corners(ind, :);
		end
		corners = cornersTmp;

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
        
        %imNewHeight = 210;
        %imNewWidth = 280;
        
		cornersNew = [         1,           1; 
					  imNewWidth,           1;
					  imNewWidth, imNewHeight;
							   1, imNewHeight];

		% Compute the homography matrix
		corners = corners';
		cornersNew = cornersNew';
        
        disp(cornersNew);
        
		h = ComputeHNorm(cornersNew, corners);
        
        % TEST
        
        %paperWidth = 210;
        %paperHeight = 280;
       
		%trimPercentage = 10;
		%eps = 1 / (1 - trimPercentage/100) - 1;
		%trimEdges = [
		%	1+eps, 0,     -eps*imNewWidth/2;
		%	0,     1+eps, -eps*imNewHeight/2;
		%	0,     0,     1;
		%];
    
        % disp(trimEdges);

        %disp(h);
        %h = trimEdges * h;
        %disp(h);
        
        % end TEST
        

		% Apply it to the original image
		tform = projective2d(h');
		imNew = imwarp(im, tform);
        
        size(imNew)
        imNew = imcrop(imNew, [600 600 3692 2498]);
        

		% Plot the results
		%subplot(2, 2, 1);
		%imshow(im); title('Original');
		%subplot(2, 2, 2);
		%imshow(foregroundMask); title('Foreground Mask');
		%subplot(2, 2, 3);
		%imshow(imGray); title('Lines & Corners');
		%hold on;
        
		%for k = 1:length(lines)
		%   xy = [lines(k).point1; lines(k).point2];
		%   plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'green');
		%   plot(xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
		%   plot(xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'red');
        %end
        
		%scatter(corners(1, :), corners(2, :), 250, 'w');
		%subplot(2, 2, 4);
		%imshow(imNew); title('Result');
        
        transformed = imNew;
end

function H2to1 = ComputeHNorm(p1, p2)
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
