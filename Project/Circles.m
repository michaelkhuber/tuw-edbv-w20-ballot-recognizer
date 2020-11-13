function ballotCircles = Circles(ballot)
% TODO - Implement this function
ballotCircles = {};
    for i = [1,2,3,4,5,6,7,8,9]
        ballotCircles{i} = imread("resources/circle_example.png");
    end
end