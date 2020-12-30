im = imread('coins.png');

[centers1, radii1] = imfindcircles(im, [15 40]);
[centers, radii] = FindCircles(im, [15 40], 15);

imshow(im);
hold on;
viscircles(centers, radii,'EdgeColor','b');
hold off
