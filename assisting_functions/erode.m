% ERODE erodes binarized image with given structuring element 
%
% Author:
%   Marie-Therese Wiedhalm
%
% Source:
%   Self
%
% Inputs:
%   im:                        The input image 
%   mask:                      Structuring element used as erosion mask  
%
% Output:
%   eImg:                      eroded image  

function [eImg] = erode(img, mask)
    eImg = ordfilt2(img, 1, mask);
end
