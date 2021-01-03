function [y] = binarize(img, thr)
    x=img;
    x=rgb2gray(x);
    [a,b]=size(x);
    y=zeros(a,b); 
    y(x>thr) = 255;  
end