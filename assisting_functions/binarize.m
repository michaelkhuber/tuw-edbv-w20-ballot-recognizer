% BINARIZE takes an image and transforms every pixel with a value bigger than the input threshold to white (255) and every other pixel to black (0) 
%
% Author:
%   Marie-Therese Wiedhalm
%
% Source:
%   Self
%
% Inputs:
%   im:                        The input image 
%   thr:                       Threshhold 
%
% Output:
%   y:                         binarized image 

function [y] = binarize(img, thr)
    x=img;
    x=rgb2gray(x);
    [a,b]=size(x);
    y=zeros(a,b); 
    y(x>thr) = 255;  
end
