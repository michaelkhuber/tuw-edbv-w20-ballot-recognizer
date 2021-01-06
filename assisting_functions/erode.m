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
%     [ai,bi]=size(img);
%     [am,bm]=size(mask);
%     deltaa = floor(am/2);
%     deltab = floor(bm/2);
%     out = img;
%     for i = ceil(am/2):ai-deltaa
%         for j = ceil(bm/2):bi-deltab
%             fenster = img(i-deltaa:i+deltaa,j-deltab:j+deltab);
%             fenster = fenster(logical(mask));
%             out(i,j) = min(fenster);
%         end
%     end
%     eImg = out;

    eImg = ordfilt2(img, 1, mask);
end
