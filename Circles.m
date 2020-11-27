function ballotCircles = Circles(ballot, ballotFilename)
    showPlot = false;
    savePlot = true;
    
    %get greyscale image
    ballot = im2gray(ballot);

    %find 10 "strongest" circles, area 15-60 maybe has to be adjusted (o in text may be recognized as circle) 
    [centers, radii, metric] = imfindcircles(ballot,[30 60]);
    
    if(showPlot)
        f = figure(2);
        clf('reset');
    elseif(savePlot)
        f = figure('visible','off');
    end
    
    if(showPlot || savePlot) 
        imshow(ballot);
        viscircles(centers, radii,'EdgeColor','b');
    end
    
    % if no circles were found, return to avoid any further errors
    if(isempty(centers))
        ballotCircles = [];
        return
    end

    %sort circles by y coordinate  
    centersSorted = sortrows(centers,2); 
    u = centersSorted(:,1);
    v = centersSorted(:,2);

    %calculate boundingboxes
    for j = 1:length(u)
        x0 = round(u(j) - radii(j));
        x1 = round(u(j) + radii(j)); 
        y0 = round(v(j) - radii(j));
        y1 = round(v(j) + radii(j));

        %create boundingboxes of circles 
        rectangle('Position', [x0 y0 radii(j)*2 radii(j)*2]);
        ballotCircles{j} = ballot(y0:y1,x0:x1);
    end
        
    if(savePlot)
        print(f,strcat("resources/results/Circles_", ballotFilename),'-dpng','-r700'); 
    end
end
