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

    %figure
    %imshow(binary_image);

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

function [area, area_size] = floodFillLabel(startX, startY, binary_image, area, area_size)
    global LABEL_MAT;
    
    nX = [1,1,1,0,0,-1,-1,-1];
    nY = [-1,0,1,-1,1,-1,0,1];
    
    % preallocate pixel stack
    stackX = zeros(size(binary_image,1) * size(binary_image,2), 1);
    stackY = zeros(size(binary_image,1) * size(binary_image,2), 1);
    
    LABEL_MAT(startX, startY) = 1;
    % set top of stack to starting pixel
    stackX(1) = startX;
    stackY(1) = startY;
    % set stack pointer
    SP = 2;
    
    while SP > 1
        area_size = area_size + 1;
        % use pixel on top of stack
        area(stackX(SP-1), stackY(SP-1)) = 1;
        
        x = stackX(SP-1) + nX;
        y = stackY(SP-1) + nY;
    
        %check if out of bounds
        isValid = (x > 0) & (y > 0) & (x <= size(binary_image, 1)) & (y <= size(binary_image,2));
        
        % remove pixel from stack
        stackX(SP) = 0;
        stackY(SP) = 0;
        SP = SP - 1;
        
        for i = find(isValid)
            if (binary_image(x(i), y(i)) == 1) && LABEL_MAT(x(i), y(i)) == 0
                LABEL_MAT(x(i), y(i)) = 1;
                % push pixel to top of stack
                stackX(SP) = x(i);
                stackY(SP) = y(i);
                SP = SP + 1;
            end
        end
    end
end
