function Main()
    original = imread("resources/ballots/4C.jpg");
    transformed = Transform2(original);    
    
    %subplot(1,2,1);
    %imshow(original); title("Original");
    %subplot(1,2,2);
    %imshow(transformed); title("Transformed");
    
end