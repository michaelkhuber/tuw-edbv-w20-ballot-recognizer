function markedCircleIndices = CheckMark(ballotCircles)
    global showPlot;
    global savePlot;
    showPlot = false;
    savePlot = false;
    markedCircleIndices = [];
   
    for k=1:length(ballotCircles)
        circle_input = ballotCircles{k};
                                  
        background = imbinarize(circle_input);
        circle = imcomplement(background);
        circle(:, 1) = 255;
        circle(:, end) = 255;
        circle(1, :) = 255;
        circle(end, :) = 255; 
        
        [~, biggest] = CountComponents(background);
        [num_components, ~] = CountComponents(circle);
        
        imshow(circle);
        
        biggest
        num_components
        
        if (num_components > 1) 
            markedCircleIndices(end+1) = k;
        elseif(biggest < 8000)
            markedCircleIndices(end+1) = k;
        end
    end
end