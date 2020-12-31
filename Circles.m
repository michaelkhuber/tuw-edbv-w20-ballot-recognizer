function ballotCircles = Circles(ballot, ballotFilename)
    global showPlot;
    global savePlot;
    showPlot = false;
    savePlot = false;
    
    global pltM;
    global pltN;
    pltM = 1;
    pltN = 1;
    pltCount = 1;
    
    if(showPlot)
        f = figure(2);
        clf('reset');
    elseif(savePlot)
        f = figure('visible','off');
    end
    
    %get greyscale image
    ballot = im2double(ballot);
    ballot = im2gray(ballot);
    
    %crop image so that only left circles are left 
    ballot=ballot(300:2150,1:500); 
    
    %find 10 "strongest" circles, area 15-60 maybe has to be adjusted (o in text may be recognized as circle) 
    % [centers, radii, metric] = imfindcircles(ballot,[30 80],'ObjectPolarity','bright','Sensitivity',0.90);
    [centers, radii] = FindCircles(ballot, [45 55]);
    [centers, radii] = SortCircles(centers, radii);
    
    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        hold on;
        imshow(ballot);
        viscircles(centers, radii,'EdgeColor','b');
        title("Hough Circles");
    end
    
    % if no circles were found, avoid any further errors
    if(isempty(centers))
        ballotCircles = [];
    else
        %find biggest circle and find circles within 0.9 % deviation
        rmax = max(radii);
        i = radii > rmax * 0.8;
        centers = centers(i, :);
        radii = radii(i);

%         centers = centers(1:10,:);
%         radii = radii(1:10,:);
        
        if(showPlot || savePlot) 
            viscircles(centers, radii,'EdgeColor','r');
            hold off;
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
            
            if(x0 > size(ballot, 2)), x0  = size(ballot, 2); end
            if(x1 > size(ballot, 2)), x1  = size(ballot, 2); end
            if(y0 > size(ballot, 1)), y0  = size(ballot, 2); end
            if(y1 > size(ballot, 1)), y1  = size(ballot, 2); end
            
            if(x0 < 1), x0 = 1; end
            if(y0 < 1), y0 = 1; end

            %create boundingboxes of circles 
            rectangle('Position', [x0 y0 radii(j)*2 radii(j)*2]);
            ballotCircles{j} = ballot(y0:y1,x0:x1);
        end
    end
        
    if(savePlot)
        print(f,strcat("resources/results/Circles_", ballotFilename),'-dpng','-r700'); 
        clf(f);
        clear f;
    end
end