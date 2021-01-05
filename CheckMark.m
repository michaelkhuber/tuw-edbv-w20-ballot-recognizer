function markedCircleIndices = CheckMark(ballotCircles)
    global showPlot;
    global savePlot;
    showPlot = true;
    savePlot = false;
    
    scan_dpixel = 5;
    
    markedCircleIndices = [];
   
    for k=1:length(ballotCircles)
        circle_input = ballotCircles{k};
                                  
        background = imbinarize(circle_input);
        circle = imcomplement(background);
        circle(:, 1) = 255;
        circle(:, end) = 255;
        circle(1, :) = 255;
        circle(end, :) = 255; 
        
        [~, ~, biggest_size] = CountComponents(background);
        [num_components, biggest, ~] = CountComponents(circle);
        
        if (num_components > 1) 
            markedCircleIndices(end+1) = k;
        elseif(biggest_size < 8000)
            markedCircleIndices(end+1) = k;
        else
           % split horizontally & count components, if there is only one 
           % empty circle, it should stay with 1 component          
           %for j=scan_dpixel: (size(circle, 1) / scan_dpixel)
           %    tmp_circle = circle;
           %    tmp_circle(2:end-1, j * scan_dpixel) = 0;
           %    [num_components, ~, ~] = CountComponents(tmp_circle);
               
           %    if (num_components > 1) 
           %       markedCircleIndices(end+1) = k;
           %       break;
           %    end
               
           %    tmp_circle = circle;
           %    tmp_circle(j * scan_dpixel, 2:end-1) = 0;
           %    [num_components, ~, ~] = CountComponents(tmp_circle);
           %    
           %    if (num_components > 1) 
           %       markedCircleIndices(end+1) = k;
           %       break;
           %    end
           %end          
        end
    end
end