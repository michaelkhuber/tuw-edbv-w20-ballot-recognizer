function [num_components, biggest, biggest_size] = CountComponents(binary_image)
% CountComponents Counts the white components in the given image via 8-way
%
% Author: 
%   Markus Hinkel
%
% Sources:
%   Haralick, Robert M., and Linda G. Shapiro, Computer and Robot Vision, Volume I, Addison-Wesley, 1992, pp. 28-48.
%
%   (https://de.mathworks.com/help/images/ref/bwconncomp.html), 
%   Official Matlab Documentation, retrieved 02.01.2020
%
% Inputs:
%   binary_image:   A binary image to count the components in
% Outputs:
%   num_components: integer with the number of found components
%   biggest:   An Integer with the pixelnumber of the biggest found
%   component
%
    num_components = 0;
    biggest = zeros(size(binary_image));
    biggest_size = 0;
    
    global LABEL_MAT;
    LABEL_MAT = zeros(size(binary_image));
    
    for i = 2  : size(binary_image, 1) -1 
       for j = 2 : size(binary_image, 2) -1
          % check if pixel is unlabeled and white, then increase component
          % number and label all adjecent pixels via floodfill
          if (LABEL_MAT(i,j) == 0) && (binary_image(i,j) == 1)
              num_components = num_components +1;
              [area, area_size] = floodFillLabel(i,j,binary_image,zeros(size(binary_image)),0);
              
              if area_size > biggest_size
                  biggest = area;
                  biggest_size = area_size;
              end
          end           
       end
    end    
end

function [res_area, res_area_size] = floodFillLabel(startX, startY, binary_image, area, area_size)
    % Mark current pixel as labeled
    global LABEL_MAT;
    LABEL_MAT(startX, startY) = 1;
    res_area_size = area_size + 1;
    res_area = area;
    res_area(startX, startY) = 1;
    
    nX = [1,1,1,0,0,-1,-1,-1];
    nY = [-1,0,1,-1,1,-1,0,1];
    
    % Check every neighbour
    for i = 1 : 8
        x = startX + nX(i);
        y = startY + nY(i);

        if (x > 0 && y > 0 && x <= size(binary_image, 1) && y <= size(binary_image,2))
            if (binary_image(x,y) == 1) && (LABEL_MAT(x,y) == 0)
                [ar, sz] = floodFillLabel(x,y, binary_image, area, area_size);
                res_area_size = res_area_size + sz;
                res_area = res_area | ar;
            end
        end
    end
end
