% SE creates an structuring element for Dilation and Erosion 
%    
%
% Author:
%   Marie-Therese Wiedhalm
%
% Source:
%   Self
%
% Inputs:
%   string:                    option if strucutring element should be a 'line' or an 'octagon' 
%   a:                         length of line or radius of octagon 
%.  b:                         degree of line 
%
% Output:
%   seOut:                     binary structuring element   

function [seOut] = se(string,a,b)
    switch string
        case 'line'
            length = a;
            deg = b;
            length = floor(length/2)*2+1;
            if deg == 0
                seOut = ones(1,(length));
            elseif deg == 90
                seOut = ones((length),1);
            end
        case 'octagon'
            r = a;
            if mod(r,3) == 0
                mask = ones(r*2+1);
                for i = 1:(2/3*r)
                    mask(i,1:2/3*r-i+1) = 0;
                    mask(i,r*4/3+2+i-1:2*r+1) = 0;
                    mask(2*r+1+1-i,1:2/3*r-i+1) = 0;
                    mask(2*r+1+1-i,r*4/3+2+i-1:2*r+1) = 0;
                end
                seOut = mask;
            else
                error('For octagon structuring elements, the SIZE input must be a nonnegative multiple of 3.') 
            end
    end
end
