function [normalizedImage, pltCount] = Normalize(image, blurStrength, saturationBased, pltCount)
% NORMALIZE processes an image such that a white foreground (e.g. white
% ballot paper) stands out in relation to the rest of the image, and afterwards
% gauss-blurs the image. 
% In this implementation, foreground-seperation is done by either blacking out 
% the background entirely, or by blurring it very strongly. In the majority of cases, a 
% histogram analysis of the saturation of the image is done. The white foreground has
% very low saturation, and its saturation distribution is thus usually seperate from
% the background saturation distribution. If the background is very dark (and thus also has
% low saturation), the brightness of the image is returned instead. If the
% background is very light and close to white (and thus again has low saturation), this
% function will most likely not yield good results, as the seperation of
% background and foreground is then quite hard with the described methods.
%
% Author:
%   Richard Binder
%
% Source:
%   Self
%
% Inputs:
%   image:              The image to process
%   blurStrength:    The sigma value for the gauss blur
%
% Output:
%   normalizedImage:    resulting normalized image

    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    
    if saturationBased
        [nobg, saturationSuccess, pltCount] = removeBackground(image, pltCount, 0.15, 0.01, false, "Saturation Histogram Analysis");
        if ( saturationSuccess )        
            nobg = im2double(nobg);
            nobg = rgb2gray(nobg);
            blurredImg = imgaussfilt(nobg, blurStrength);
            normalizedImage = blurredImg;
            return;
        end

        [nobg, saturationSuccess, pltCount] = removeBackground(image, pltCount, 0.15, 0.3, true, "Extended Saturation Histogram Analysis");
        if ( saturationSuccess )
            nobg = im2double(nobg);
            nobg = rgb2gray(nobg);    

            chroma = getChroma(image);
            white = sqrt(3*double(255)^2);
            chroma = sqrt(double(chroma(:,:,1)).^2 + double(chroma(:,:,2)).^2 + double(chroma(:,:,3)).^2) ./ white;
            chroma = adapthisteq(chroma,'clipLimit',0.02,'Distribution','rayleigh');

            chroma(nobg < 0.1) = chroma(nobg < 0.1) - 3 * std(chroma(:));
            chroma( chroma < 0.01 ) = 0.0;
            gray = rgb2gray(im2double(image));
            %gray = sqrt( chroma .* gray );
            gray(nobg < 0.1) = gray(nobg < 0.1) - 2 * std(gray(:));
            blurredImg = imgaussfilt(gray, 100);
            blurredImg(nobg > 0.1) = gray(nobg > 0.1);
            blurredImg = imgaussfilt(blurredImg, blurStrength/2);

            if(showPlot || savePlot) 
                subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
                imshow(blurredImg); title("Background Based Chromaticity Blurring");
            end
            normalizedImage = blurredImg;
            return;
        end
    end
    
    grayImg = rgb2gray(im2double(image));
    blurredImg = imgaussfilt(grayImg, blurStrength);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        imshow(blurredImg); title("Fallback To Denoised Brightness");
    end
    
    normalizedImage = blurredImg;
end

% Removes Background based on the histogram of the saturation of the input image.
% On a local minimum of the saturation histogram, a threshold is applied onto
% the image. Everything below that threshold is blacked out. If the number
% of blacked out pixels is below 20% or above 80% or if no significant
% local minimum could be found, it is assumed that the
% seperation between saturation distributions failed, and success is
% returned as false, otherwise true. maxProminence and minProminence are
% parameters for finding a local minimum and maximum that define the
% saturation distributions. For more info see description of function
% "getLocals".
function [thresholdImage, success, pltCount] = removeBackground(image, pltCount, maxProminence, minProminence, doZeroApprox, histogramTitle)
    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    
    hsv = rgb2hsv(image);    
    mask = 1.0 - hsv(:,:,2);
    
    [counts, edges] = histcounts(mask, 100);
    [locMin, locMax] = getLocals(counts, maxProminence, minProminence);
    edges = edges(1:(length(edges)-1));
    
    if( isempty(locMin) ) 
        thresholdImage = image;
        success = 0;
        return;
    end    
    
    localMinDist = edges(locMin);
    
    if( doZeroApprox )
        zeroApprox = locMin - 5 * counts(locMin) * (locMax - locMin) / (counts(locMax) - counts(locMin));
        zeroApprox = round( zeroApprox );
        if( zeroApprox < 1 )
            zeroApprox = 1;
        end
        localMinDist = edges(round(zeroApprox));
    end
    
    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        hold on;
        plot(1.0-edges,counts,'Color', 'red');
        plot(1.0-edges(locMax),counts(locMax),'^','Color', 'blue');
        plot(1.0-edges(locMin),counts(locMin),'v','Color', 'blue');
        
        title(histogramTitle);
        hold off;
    end
    
    binary = mask <= localMinDist;
    pixelsRemoved = sum(binary(:));
    
    thresholdImage = image;    
    for i = 1:3
        channel = image(:,:,i);
        channel(binary) = 0.0;
        thresholdImage(:,:,i) = channel;
    end
    
    % assume that number of removed pixels is an indicator whether the
    % background removal succeeded or not
    minRemoved = 0.2 * length(binary(:));
    maxRemoved = 0.8 * length(binary(:));
    
    if((pixelsRemoved < minRemoved) || (pixelsRemoved > maxRemoved) )
        success = 0;
    else
        success = 1;
        if(showPlot || savePlot) 
            subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
            imshow(thresholdImage); title("Saturation Based Background Removal");
        end
    end
end

% computes a local minimum and local maximum of a saturation histogram "counts".
% There is several conditions that this local maximum and minimum need to fulfill:
%   - The local  maximum needs to be at least maxProminence times the
%   global maximum of counts.
%   - The local minimum needs to be less than minProminence times that
%   local maximum. This ensures that this local minimum is significantly
%   smaller than the local maximum and indicates the seperation point between
%   two different saturation distributions.
%   - The local minimum needs to be to the right of the local maximum. This
%   ensures that the local maximum is the peak of the left most saturation
%   distribution and therefore indicates a white (low saturated) foreground.
function [locMin, locMax] = getLocals(counts, maxProminence, minProminence)
    maxProminence = maxProminence * max(counts(:));

    %get maxima
    maxima = [];
    for i = 2:(length(counts)-1)
        if( (counts(i) > counts(i-1)) && (counts(i) > counts(i+1)) )
            maxima = [maxima, i];
        end
    end
    
    %edge cases
    if( counts(1) > counts(2) )
        maxima = [1, maxima];
    end
    if( counts(length(counts)) > counts(length(counts)-1) )
        maxima = [maxima, length(counts)];
    end
    
    minima = [];
    for i = 2:(length(counts)-1)
        if( (counts(i) < counts(i-1)) && (counts(i) < counts(i+1)) )
            minima = [minima, i];
        end
    end
    
    %edge cases
    if( counts(1) < counts(2) )
        minima = [1, minima];
    end
    if( counts(length(counts)) < counts(length(counts)-1) )
        minima = [minima, length(counts)];
    end
    
    prominentMaxima = [];
    prominentMinima = [];
    for i = 1:length(maxima)
        if( counts(maxima(i)) > maxProminence )
            for j = 1:length(minima)
                if( minima(j) < maxima(i) && counts(minima(j)) < minProminence * counts(maxima(i)) )
                    prominentMaxima = [prominentMaxima, maxima(i)];
                    prominentMinima = [prominentMinima, minima(j)];
                end
            end
        end
    end
    
    %get right most prominent maximum
    locMax = max(prominentMaxima);
    
    %and its right most assigned prominent minimum
    locMin = max(prominentMinima(prominentMaxima==locMax));
end

%gets Chromaticity of a (uint8) rgb image
function chroma = getChroma(image)
    image = im2double(image);
    brightness = rgb2gray(image);
    brightness(brightness < 0.01) = 0.01;
    r = image(:,:,1) ./ brightness;
    g = image(:,:,2) ./ brightness;
    b = image(:,:,3) ./ brightness;
    chroma = cat(3, r, g, b);
    chroma = im2uint8(chroma);
end