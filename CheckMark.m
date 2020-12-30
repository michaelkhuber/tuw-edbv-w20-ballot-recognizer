function markedCircleIndices = CheckMark(ballotCircles)
    global showPlot;
    global savePlot;
    showPlot = false;
    savePlot = false;
    markedCircleIndices = [];
   
    for k=1:length(ballotCircles)
        circle_input = ballotCircles{k};
                                  
        circle = imbinarize(circle_input);
        circle = imcomplement(circle);
        circle(:, 1) = 255;
        circle(:, end) = 255;
        circle(1, :) = 255;
        circle(end, :) = 255; 
        
        CC = bwconncomp(circle);
        
        if(showPlot || savePlot)
            
        end
        
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idx] = max(numPixels);
        
        if (CC.NumObjects > 1) 
            markedCircleIndices(end+1) = k;
        elseif(biggest > 2000)
            markedCircleIndices(end+1) = k;
        end
    end
end