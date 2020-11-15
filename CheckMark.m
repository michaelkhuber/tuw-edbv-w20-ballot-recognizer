function markedCircleIndices = CheckMark(ballotCircles)
    markedCircleIndices = [];
   
    for k=1:length(ballotCircles)
        circle_rgb = ballotCircles{k};
        
        % convert to black and white
        circle = im2gray(circle_rgb);

        % convert to binary
        circle = imbinarize(circle);
        
        % hough transform to find lines
        [H,theta,rho] = hough(circle);

        P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));

        x = theta(P(:,2));
        y = rho(P(:,1));
        plot(x,y,'s','color','black');
        
        % lines have to have particular length
        lines = houghlines(circle,theta,rho,P,'FillGap',2.5,'MinLength',16);
    
        % if at least one line is on/in the circle, it was checked
        if(length(lines) > 1)
           markedCircleIndices(end+1) = k;
           
        end
        
        figure, imshow(circle_rgb), hold on
        
        for k = 1:length(lines)
            xy = [lines(k).point1; lines(k).point2];
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

        end
    end
    
end