% GAUSSFILT performs a gaussian blur on the input image with given sigma 
%
% Author:
%   Marie-Therese Wiedhalm
%
% Source:
%   Self
%
% Inputs:
%   im:                        The input image 
%   sigma:                     standard deviation   
%
% Output:
%   gImg:                      Image with gaussian blur   

function [gImg] = gaussfilt(img, sigma)
    delta = round(3*sigma);
    kernelWidth = 2*delta+1; 
    g = zeros(kernelWidth,kernelWidth); 
    
    for x=-delta:delta
      for y=-delta:delta
        g(x+delta+1,y+delta+1) = 1 / (2*pi*sigma^2) / exp((x^2+y^2) / (2 * sigma^2)); 
      end 
    end
    sumG = sum(g,'all');
    g = g/sumG;
    gImg = imfilter(img,g,'replicate'); 
end 
