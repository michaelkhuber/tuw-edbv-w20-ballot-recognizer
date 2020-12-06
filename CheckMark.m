function markedCircleIndices = CheckMark(ballotCircles)
    global showPlot;
    global savePlot;
    showPlot = false;
    savePlot = false;
    markedCircleIndices = [];
   
    for k=1:length(ballotCircles)
        circle_input = ballotCircles{k};
                         
        % convert to binary
        %circle = imbinarize(circle_input);
         
        circle = imbinarize(circle_input);
        circle = edge(circle,'canny');
                
        % hough transform to find lines
        [H,theta,rho] = hough(circle);

        P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));

        x = theta(P(:,2));
        y = rho(P(:,1));
        if(showPlot || savePlot)
            plot(x,y,'s','color','black');
        end
        
        % lines have to have particular length
        lines = houghlines(circle,theta,rho,P,'FillGap',8.0,'MinLength',45);
    
        % if at least one line is on/in the circle, it was checked
        if(length(lines) > 1)
           markedCircleIndices(end+1) = k;
           
        end
        
        %figure, imshow(circle_input), hold on
        
        %for n = 1:length(lines)
            %xy = [lines(n).point1; lines(n).point2];
            %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

        %end
    end
    
end