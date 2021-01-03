function [resizedImg] = resize(img, newSize)    
    [rows, columns, colorChannels] = size(img);
    if colorChannels == 1 
        [rows, columns] = size(img);
        [X,Y] = meshgrid(1:columns,1:rows);
        xi = linspace(1,columns,newSize(2));
        yi = linspace(1,rows,newSize(1));
        [Xi,Yi] = meshgrid(xi,yi); 
        resizedImg = uint8(interp2(X,Y,double(img),Xi,Yi));
    else    
        [X,Y,Z] = meshgrid(1:columns,1:rows,1:colorChannels);
        xi = linspace(1,columns,newSize(2));
        yi = linspace(1,rows,newSize(1));
        zi = linspace(1,colorChannels,3); 
        [Xi,Yi,Zi] = meshgrid(xi,yi,zi); 
        resizedImg = uint8(interp3(X,Y,Z,double(img),Xi,Yi,Zi));
    end 
end
