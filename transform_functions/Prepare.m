function preparedImage = Prepare(image, blurStrength, saturationBased)
% PREPARE processes an image such that a white foreground (e.g. white
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
%   preparedImage:    resulting prepared image

    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    global pltCount;
    
    if saturationBased
        [nobg, saturationSuccess] = removeBackground(image, 0.15, 0.01, false, "Saturation Histogram Analysis");
        if ( saturationSuccess )        
            nobg = im2double(nobg);
            nobg = toGray(nobg);
            blurredImg = gaussfilt(nobg, blurStrength);
            preparedImage = blurredImg;
            return;
        end

        [nobg, saturationSuccess] = removeBackground(image, 0.15, 0.05, true, "Extended Saturation Histogram Analysis");
        if ( saturationSuccess )
            nobg = im2double(nobg);
            nobg = toGray(nobg);    
            blurredImg = gaussfilt(nobg, blurStrength);
            preparedImage = blurredImg;
            return;
        end
    end
    
    grayImg = toGray(im2double(image));
    %grayImg = adapthisteq(grayImg,'NumTiles',[8 8],'ClipLimit',0.005);
    blurredImg = gaussfilt(grayImg, blurStrength);

    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        imshow(blurredImg); title("Fallback To Denoised Brightness");
    end
    
    preparedImage = blurredImg;
end

function [thresholdImage, success] = removeBackground(image, maxProminence, minProminence, doZeroApprox, histogramTitle)
% REMOVEBACKGROUND removes the background of an image based on the 
% histogram of the saturation of the input image.
% On a local minimum of the saturation histogram, a threshold is applied onto
% the image. Everything above that threshold is blacked out. If the number
% of blacked out pixels is below 20% or above 80% or if no significant
% local minimum could be found, it is assumed that the seperation between 
% saturation distributions (and thus also the background removal) failed.
%
% Author:
%   Richard Binder
%
% Source:
%   Self
%
% Inputs:
%   maxProminence:    prominence parameter for finding local maximum, see 
% description of function "getLocals".
%   minProminence:     prominence parameter for finding local minimum, see 
% description of function "getLocals".

    global showPlot;
    global savePlot;
    global pltM;
    global pltN;
    global pltCount;
    
    hsv = rgb2hsv(image);    
    mask = hsv(:,:,2);
    
    [counts, edges] = histcounts(mask, 100);
    [locMin, locMax] = getLocals(counts, maxProminence, minProminence);
    edges = edges(1:(length(edges)-1));
    
    if( isempty(locMin) ) 
        thresholdImage = image;
        success = 0;
        return;
    end    
    
    localMinDist = edges(locMin);
    
    if doZeroApprox 
        zeroApprox = locMin - 2 * (locMax - locMin) * counts(locMin) / (counts(locMax) - counts(locMin));
        zeroApprox = round( zeroApprox );
        if( zeroApprox < 1 )
            zeroApprox = 1;
        end
        localMinDist = edges(round(zeroApprox));
    end
    
    if(showPlot || savePlot) 
        subplot(pltM, pltN, pltCount); pltCount = pltCount + 1;
        hold on;
        plot(edges,counts,'Color', 'red');
        plot(edges(locMax),counts(locMax),'^','Color', 'blue');
        plot(edges(locMin),counts(locMin),'v','Color', 'blue');
        if doZeroApprox
            plot(edges(zeroApprox), counts(zeroApprox), 'o', 'Color', 'green');
        end
        
        title(histogramTitle);
        hold off;
    end
    
    binary = mask > localMinDist;
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

function [locMin, locMax] = getLocals(counts, maxProminence, minProminence)
% GETLOCALS computes a local minimum and a local maximum of a saturation histogram "counts".
% There is several conditions that this local maximum and minimum need to fulfill:
%   - The local maximum needs to be at least maxProminence times the
%   global maximum of counts.
%   - The local minimum needs to be less than minProminence times the
%   local maximum. This ensures that this local minimum is significantly
%   smaller than the local maximum and indicates the seperation point between
%   two different saturation distributions.
%   - The local maximum needs to be to the left of the local minimum. This
%   ensures that the local maximum is the peak of the left most saturation
%   distribution and therefore indicates a white (low saturated) foreground.
%
% Author:
%   Richard Binder
%
% Source:
%   Self
%
% Inputs: 
%   counts:                     an array containing histogram values
%   maxProminence:     prominence parameter for finding local maximum
%   minProminence:      prominence parameter for finding local minimum
%
% Outputs:
%   locMin:                     the index of the local minimum in the
%   counts array
%   locMax:                     the index of the local maximum in the
%   counts array


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
                if( minima(j) > maxima(i) && counts(minima(j)) < minProminence * counts(maxima(i)) )
                    prominentMaxima = [prominentMaxima, maxima(i)];
                    prominentMinima = [prominentMinima, minima(j)];
                end
            end
        end
    end
    
    %get left most prominent maximum
    locMax = min(prominentMaxima);
    
    %and its left most prominent minimum
    locMin = min(prominentMinima(prominentMaxima==locMax));
end