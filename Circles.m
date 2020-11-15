function ballotCircles = Circles(ballot)

    %find 10 "strongest" circles, area 15-60 maybe has to be adjusted (o in text may be recognized as circle) 
    imshow(ballot);
    [centers, radii, metric] = imfindcircles(ballot,[30 60]);
    viscircles(centers, radii,'EdgeColor','b');

    %sort circles by y coordinate  
    centersSorted = sortrows(centers,2); 
    u = centersSorted(:,1);
    v = centersSorted(:,2);
    length(u)

    %calculate boundingboxes
    for j = 1:length(u)
        x0 = u(j) - radii(j);
        x1 = u(j) + radii(j); 
        y0 = v(j) - radii(j);
        y1 = v(j) + radii(j);

        %create boundingboxes of circles 
        rectangle('Position', [x0 y0 radii(j)*2 radii(j)*2]);
        ballotCircles{j} = ballot(y0:y1,x0:x1);
    end
    
end
