% TOGRAY takes an RGB image and returns a grayscale image, accepts grayscale images as inputs and returns them unmodified  
%
% Author:
%   Marie-Therese Wiedhalm
%
% Source:
%   Self
%
% Inputs:
%   im:                        The input image   
%
% Output:
%   gray:                      grayscale image  

function [gray] = toGray(img)
    [a,b,c]=size(img);
    if c==3 
        gray = 0.2989 * img(:,:,1) + 0.5870 * img(:,:,2) + 0.1140 * img(:,:,3);
    else
        gray=A;
    end
end
