% RESIZE resizes an input image to new dimensions using interp3() 
%    
%
% Author:
%   Marie-Therese Wiedhalm
%
% Source:
%   Implementation is based on Image Analyst(2020). Shrink image without using imresize 
%   https://uk.mathworks.com/matlabcentral/answers/554017-shrink-image-without-using-imresize
%   MATLAB Central Answers. Retrieved December 25, 2020.

% Inputs:
%   img:                       input image 
%   newSize:                   vector of new size of image  
%
% Output:
%   resizedImg:                resized image   

function [resizedImg] = resize(img, newSize)    
    [rows, columns, colorChannels] = size(img);
    % resizing if image is grayscale 
    if colorChannels == 1 
        [rows, columns] = size(img);
        [X,Y] = meshgrid(1:columns,1:rows);
        xi = linspace(1,columns,newSize(2));
        yi = linspace(1,rows,newSize(1));
        [Xi,Yi] = meshgrid(xi,yi); 
        resizedImg = uint8(interp2(X,Y,double(img),Xi,Yi));
    % resizing if image is RGB 
    else    
        [X,Y,Z] = meshgrid(1:columns,1:rows,1:colorChannels);
        xi = linspace(1,columns,newSize(2));
        yi = linspace(1,rows,newSize(1));
        zi = linspace(1,colorChannels,3); 
        [Xi,Yi,Zi] = meshgrid(xi,yi,zi); 
        resizedImg = uint8(interp3(X,Y,Z,double(img),Xi,Yi,Zi));
    end 
end
