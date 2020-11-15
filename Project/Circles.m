function ballotCircles = Circles(ballot)
% TODO - Implement this function
%ballotCircles = {};
%   for i = [1,2,3,4,5,6,7,8,9]
%      ballotCircles{i} = imread("resources/circle_example.png");
% end

%I = rgb2gray(ballot);
%imshow(I);

%find 10 "strongest" circles, area 15-60 maybe has to be adjusted (o in text may be recognized as circle) 
[centers, radii, metric] = imfindcircles(I,[15 60]);
viscircles(centers, radii,'EdgeColor','b');

%sort circles by y coordinate  
centersSorted = sortrows(centers,2); 
u = centersSorted(:,1);
v = centersSorted(:,2);
length(u)

%calculate boundingboxes
for j = 1:length(u)
x0 = u(j) - radii(j);
x1 = u(j) + radii(j); 
y0 = v(j) - radii(j);
y1 = v(j) + radii(j);

%create boundingboxes of circles 
rectangle('Position', [x0 y0 radii(j)*2 radii(j)*2]);
    switch j
        case 1
            P1 = A(y0:y1,x0:x1);
        case 2
            P2 = A(y0:y1,x0:x1);
        case 3
            P3 = A(y0:y1,x0:x1);
        case 4
            P4 = A(y0:y1,x0:x1);
        case 5
            P5 = A(y0:y1,x0:x1);
        case 6
            P6 = A(y0:y1,x0:x1);
        case 7
            P7 = A(y0:y1,x0:x1);
        case 8
            P8 = A(y0:y1,x0:x1);
        case 9
            P5 = A(y0:y1,x0:x1);
        case 10
            P10 = A(y0:y1,x0:x1);
    end
end
    
end
