function [lines, lines2, pltCount] = HoughTransform(maskedImage, pltCount)
        global showPlot;
        global savePlot;
        
        maxPeaks = 50;
        rhoPrecision = 5;
        
        maxRho = sqrt(size(maskedImage,1)^2 + size(maskedImage,2)^2);
        maxRho = ceil(maxRho);
        rhoInterval = [-maxRho : rhoPrecision : maxRho];
        
		% Use Hough transform to find vertical lines
        thetaInterval = [-45.0 : 0.25 : 44.75];
		%[H, theta, rho] = hough(maskedImage, 'RhoResolution', rhoPrecision, 'Theta', thetaInterval);
		[H, theta, rho] = Hough(maskedImage, rhoInterval, thetaInterval);
        P = Houghpeaks(H, maxPeaks);
		%lines = houghlines(maskedImage, theta, rho,P, 'MinLength', 10);
		lines = Houghlines(maskedImage, thetaInterval, rhoInterval, P);
        
		% Use Hough transform to find horizontal lines
        thetaInterval2 = [-90.0 : 0.25 : -44.75];
        thetaInterval2 = [thetaInterval2, 45.0 : 0.25 : 89.75];
		%[H2, theta2, rho2] = hough(maskedImage, 'RhoResolution', rhoPrecision, 'Theta', thetaInterval);
		[H2, theta2, rho2] = Hough(maskedImage, rhoInterval, thetaInterval2);
        P2 = Houghpeaks(H2, maxPeaks);
		%lines2 = houghlines(maskedImage, theta2, rho2,P2, 'FillGap', 500, 'MinLength', 500);
		lines2 = Houghlines(maskedImage, thetaInterval2, rhoInterval, P2);
        
        if(showPlot || savePlot)
             pltCount = plotHoughTransform(H, theta, rho, P, H2, theta2, rho2, P2, pltCount);
        end
end





function [ H, thetaInterval, rhoInterval ] = Hough( maskedImage, rhoInterval, thetaInterval )
    H = zeros(length(rhoInterval),length(thetaInterval));
    [y, x] = find(maskedImage);

    for thetaIndex = 1:length(thetaInterval)
        theta = thetaInterval(thetaIndex);
        rho = x*cosd(theta) + y*sind(theta);
        rhoIndices = getRhoIndex(rho, rhoInterval);
        thetaIndices =repmat(thetaIndex, [length(rhoIndices) 1]);
        H = H + accumarray([rhoIndices, thetaIndices],1,size(H));
    end
end

function rhoIndex = getRhoIndex(rho, rhoInterval)
    maxRho = max(rhoInterval);
    minRho = min(rhoInterval);
    rhoIndex = (rho + abs(minRho) ) / ( abs(minRho) + maxRho ) * length(rhoInterval);
    rhoIndex = floor(rhoIndex) + 1;
end



function P = Houghpeaks(H, maxPeaks)
    peaks = ordfilt2(H, 9, ones(3, 3));
    peaksBinary = (H == peaks) & (H > 0);
    
    threshold = 0.5 * max(H(:));
    peaksBinary = peaksBinary & (peaks > threshold);
    
    [rho, theta] = find(peaksBinary);
    P = [rho, theta];
end


function lines = Houghlines(maskedImage, thetaInterval, rhoInterval, P)
    lines = [];
    for i = [1 : size(P, 1)]
        peak = P(i,:);
        line.rho = rhoInterval(peak(1));
        line.theta = thetaInterval(peak(2));
        
        [line.point1, line.point2] = findEndpoints(maskedImage, line, rhoInterval);
        
        lines = [lines, line];
    end
end


function [point1, point2] = findEndpoints(maskedImage, line, rhoInterval)
    [y, x] = find(maskedImage);
    
    rho = x*cosd(line.theta) + y*sind(line.theta);
    rhoIndices = getRhoIndex(rho, rhoInterval);
    lineRhoIndex = getRhoIndex(line.rho, rhoInterval);
    
    coords = [x, y];
    linePixels = coords(rhoIndices == lineRhoIndex, :);
    linePixels = sortrows(linePixels, [1 2]);
    
    point1 = linePixels(1,:);
    point2 = linePixels(length(linePixels),:);
end


function pltCount = plotHoughTransform(H, theta, rho, P, H2, theta2, rho2, P2, pltCount)
    global pltM;
    global pltN;
    
    %plot Vertical Hough transformation
    subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
    imshow(imadjust(rescale(H)),'XData',theta,'YData',rho, 'InitialMagnification','fit');
     axis on, axis normal, hold on;
     xlabel('\theta'), ylabel('\rho');
    plot(theta(P(:,2)),rho(P(:,1)),'s','color','white'); %plot peaks
    colormap(gca,hot);
     title('Hough Transform Vertical');
    hold off;

    %plot Horizontal Hough transformation
    subplot(pltM, pltN, pltCount);  pltCount = pltCount + 1;
    imshow(imadjust(rescale(H2)),'XData',theta2,'YData',rho2);
     axis on, axis normal, hold on;
     xlabel('\theta'), ylabel('\rho');
    plot(theta2(P2(:,2)),rho2(P2(:,1)),'s','color','white'); %plot peaks
    colormap(gca,hot);
     title('Hough Transform Horizontal');
    hold off;
end