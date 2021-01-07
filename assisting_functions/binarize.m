% BINARIZE takes an image and marks every pixel with a value bigger than graythresh(img) 
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

function y = binarize(img)
    y = img > graythresh(img);
end
