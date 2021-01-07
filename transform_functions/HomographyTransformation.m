function imNew = HomographyTransformation(image, corners, step)
% HOMOGRAPHYTRANSFORMATION determines the projection of an image from the
% 3d world to the 2d plane assuming that corners are 4 points on the 2d
% plane that would make up a perfect rectangle in the 3d world. The image
% is then transformed by undoing the projection such that the corners make
% up a perfect rectangle again.
%
% Author:
%   Jakob, Richard
%   The code is partially taken from
%   https://github.com/jasonfly07/matlab_ws/blob/master/document_scanner/main.m
%
% Source:
%   The function is based on solving the linear equation for Homography
%   Estimation in which 3d points are projected onto the 2d plane by
%   multiplying them with a projection matrix. More on the topic can be
%   found at:
%   http://cseweb.ucsd.edu/classes/wi07/cse252a/homography_estimation/homography_estimation.pdf
%
% Inputs:
%   image:      the image to transform
%   corners:   4 points that would make up a perfect rectangle in the 3d
%   world
%
% Output:
%   imNew:      the transformed image

    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    global pltCount;

    % Measure the (skewed) widths & heights of the quadrangle
    heightL = norm(corners(1,:) - corners(4,:));
    heightR = norm(corners(2,:) - corners(3,:));
    widthT = norm(corners(1,:) - corners(2,:));
    widthB = norm(corners(3,:) - corners(4,:));
    
    % cancel if points are too close to each other (results in Matlab
    % crashing when performing imwarp)
    if (heightL < 300 || heightR < 300 || widthT < 300 || widthB < 300)
        error("Points of quadrangle too close to each other, cannot perform imwarp");
    end
    
    % We know the dimensions of our ballot paper and ballot table, so we
    % can just use them to define the dimensions of the new rectangle
    if step == 1
        ballotHeight = 2100;
        ballotWidth = 2970;
    elseif step == 2
        ballotHeight = 2150;
        ballotWidth = 3520;
    end
    
    % compute the skewfactor, which should give the ratio between
    % horizontal and vertical points
    skewFactor = ( widthT * widthB ) / ( heightL * heightR );
    
    %depending on the skewFactor, we decide if the horizontal or vertical
    %lines are the longer ones.
    if skewFactor > 1.1
        imNewWidth = ballotWidth;
        imNewHeight = ballotHeight;
    else
        imNewWidth = ballotHeight;
        imNewHeight = ballotWidth;
    end

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
    [imNew, RB] = imwarp(image, tform);

    if(showPlot || savePlot) 
        % Plot the results
        pltCount = pltCount + 1; subplot(pltM, pltN, pltCount);
        t = 'Transformed';
        imshow(imNew); title([num2str(pltCount), '. ', t]);
    end

    %Crop the image to to the 4 corners of the ballot
    upper_left  = corners([1,2],1)';
    lower_right = corners([1,2],3)';

    [lower_right_X,lower_right_Y] = transformPointsForward(tform, lower_right(1), lower_right(2));
    [X1, Y1] = worldToIntrinsic(RB, lower_right_X, lower_right_Y);

    [upper_left_X, upper_left_Y] = transformPointsForward(tform, upper_left(1), upper_left(2));
    [X2, Y2] = worldToIntrinsic(RB, upper_left_X, upper_left_Y);

    imNew = imcrop(imNew, [X2 Y2 X1-X2 Y1-Y2]);
end

function H2to1 = ComputeHNorm(p1, p2)
% COMPUTEHNORM computes the homography projection matrix between two point sets
%
% Author:
%   Jakob
%   The code is largely taken from
%   https://github.com/jasonfly07/matlab_ws/blob/master/document_scanner/ComputeHNorm.m
%
% Inputs:
%   p1:      the first point set, made up of 4 points
%   p2:      the second point set, made up of 4 points
%
% Output:
%   H2to1:  the projection matrix

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
% NORM_MATRIX computes the norm matrix of a point set
%
% Author:
%   Jakob
%   The code is largely taken from
%   https://github.com/jasonfly07/matlab_ws/blob/master/document_scanner/ComputeHNorm.m
%
% Inputs:
%   p2:      the point set
%
% Output:
%   Tnorm:  the norm matrix

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