% DILATE dilates binarized image with given structuring element 
%
% Author:
%   Marie-Therese Wiedhalm
%
% Source:
%   Self
%
% Inputs:
%   im:                        The input image 
%   mask:                      Structuring element used as dilation mask  
%
% Output:
%   dImg:                      dilateted image  

function [dImg] = dilate(img, mask)
    dImg = ordfilt2(img, length(find(mask)), mask);
end
