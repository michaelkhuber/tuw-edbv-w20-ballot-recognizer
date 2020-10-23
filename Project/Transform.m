function transformedBallot = Transform(ballot, template)
    original = template;
    distorted = ballot;

    distorted = rgb2gray(distorted);
    original = rgb2gray(original);

    % Detect and extract features from the original and the transformed images.
    ptsOriginal  = detectSURFFeatures(original);
    ptsDistorted = detectSURFFeatures(distorted);
    [featuresOriginal,validPtsOriginal] = extractFeatures(original,ptsOriginal);
    [featuresDistorted,validPtsDistorted] = extractFeatures(distorted,ptsDistorted);
    
    % Match and display features between the images.
    index_pairs = matchFeatures(featuresOriginal,featuresDistorted);
    matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
    matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));
    figure 
    showMatchedFeatures(original,distorted,matchedPtsOriginal,matchedPtsDistorted)
    title('Matched SURF Points With Outliers');
    
    % Exclude the outliers, estimate the transformation matrix, and display the results.
    [tform,inlierIdx] = estimateGeometricTransform2D(matchedPtsDistorted,matchedPtsOriginal,'similarity');
    inlierPtsDistorted = matchedPtsDistorted(inlierIdx,:);
    inlierPtsOriginal  = matchedPtsOriginal(inlierIdx,:);

    figure 
    showMatchedFeatures(original,distorted,inlierPtsOriginal,inlierPtsDistorted)
    title('Matched Inlier Points')
    
    % Use the estimated transformation to recover and display the original image from the distorted image.
    outputView = imref2d(size(original));
    recovered = imwarp(distorted,tform,'OutputView',outputView);
    figure 
    imshow(recovered); 
    title('Recovered Image');

    transformedBallot = recovered;
end