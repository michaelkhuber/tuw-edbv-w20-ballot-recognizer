function Main()
    I = imread("resources/ballots/4C.jpg");
    I2 = Transform(I);    
    
    subplot(1,2,1);
    imshow(I); title("Original");
    subplot(1,2,2);
    imshow(I2); title("Transformed");
end