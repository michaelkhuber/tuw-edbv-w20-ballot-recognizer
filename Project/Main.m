% - Our pipeline can return 4 different outputs:
% - ballot is invalid
% - ballot's validity is unidentifiable
% - ballot is inaproppriately marked
% - ballot is valid and apropproiately marked -> return marked choice

% ###########
% ### MAIN ###
% ###########
function choices = Main()
    % - Read in a ballot template image and possible template choices
    [template, templateChoices] = Templ();

    % - Read in all Ballot Filenames from the Ballot Folder
    ballotFilenames = BallotFilenames();

    % - Preallocate choices array
    choices = strings(length(ballotFilenames));
    % - Do something for every ballot image
    i = 1;
    for ballotFilename = ballotFilenames
        % Figure out the marked choice for a ballot image
        choice = Pipeline(template, templateChoices, ballotFilename);
        choices(i) = choice;
    end
end

% ###############
% ### PIPELINE ###
% ###############
function choice = Pipeline(template, templateChoices, ballotFilename)
        % - For each Ballot, go through the following steps:
        
        % - STEP 1
        % - Read in ballot image
        % - Prepare the image for the template matching
        ballotImg = Read(ballotFilename);
        imshow(ballotImg);
        
        % - STEP 2
        % - Match template to ballot
        % - Includes any necessary gemeotric tansformations
        % - If a match cannot be made, declare the ballot invalid -> cancel the
        % pipeline and return
        % - this step can maybe be split into two functions instead of just one?
        transformedBallot = Transform(ballotImg, template);

        % - STEP 3
        % - Match circles in ballot
        % - Associate each circle with its appropriate choice by
        % outputting circles in correct order
        % - If the right amount of circles cannot be found, declare the ballot's
        % validity as unidentifiable -> cancel the pipeline and return
        ballotCircles = Circles(transformedBallot);
        
        % - STEP 4
        % - Figure out which circle is marked
        % - If a marked circle is found, declare the ballot valid and return 
        % - If no circle is marked, declare the ballot inaproppriately marked ->
        % cancel the pipeline and return
        markedCircleIndex = CheckMark(ballotCircles);
        choice = templateChoices(markedCircleIndex);
end

