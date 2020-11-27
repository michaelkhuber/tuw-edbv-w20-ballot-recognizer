% - Our pipeline can return 4 different outputs:
% - ballot is invalid
% - ballot's validity is unidentifiable
% - ballot is inaproppriately marked
% - ballot is valid and apropproiately marked -> return marked choice


%% ##########
%  ## MAIN ##
%  ##########
function ballotTable = Main()
    % - Read in ordered ballot template choices
    templateChoices = Templ();

    % - Read in all Ballot Filenames from the Ballot Folder
    ballotFilenames = BallotFilenames();
    % - Manually choose Filenames (meant for debugging)
    ballotFilenames = ["resources/ballots/2A.jpg", "2A.jpg"];
    
    numBallots = size(ballotFilenames,1);
    % - Preallocate validity array
    % - Each entry can have take one of the three values:
    % - valid
    % - invalid
    % - unidentified
    validity = strings([numBallots, 1]);
    % - Preallocate choices array
    % - Each entry contains the marked choices in a ballot (if any)
    choices = strings([numBallots, 1]);
    % - Preallocate errors array
    % - Each entry contains the reason why the ballot couldnt be validated
    errors = strings([numBallots, 1]);
    
    % - Do something for every ballot image
    for i = 1:numBallots
        % Figure out the marked choice for a ballot image (if it is valid)
        [validity(i), choices(i), errors(i)] = Pipeline(templateChoices, ballotFilenames(i,:));
    end
    
    index = (1:numBallots)';
    colnames = ["index", "filename", "validity", "choices", "errors"];
    names = ballotFilenames(:,2);
    ballotTable = table(index, names, validity, choices, errors, 'VariableNames', colnames);
    writetable(ballotTable, "result.csv")
end

%% ##############
%  ## PIPELINE ##
%  ##############
function [validity, choice, error] = Pipeline(templateChoices, ballotFilename)
    error = "";
    try
        % - Implemented as suggested in the file
        % "Konzept_Wahlzettel_Erkennung.pdf" -> Point 5: Methodik
        % - For each Ballot, go through the following steps:
        
        %% - STEP 1
        %  - Read in ballot image
        ballotImg = Read(ballotFilename(1));
        
        %  - Prepare the image for circle matching
        % preparedBallot = Prepare(ballotImg);
        
        %% - STEP 2
        %  - Normalize and transform the image
        transformedBallot = Transform(ballotImg, ballotFilename(2));
        
        %% - STEP 3
        %  - Match circles in ballot
        %  - output all circles in correct order
        ballotCircles = Circles(transformedBallot, ballotFilename(2));
        
        %  - If the right amount of circles cannot be found, declare the ballot's
        %  validity as unidentifiable -> cancel the pipeline and return
        if length(ballotCircles) ~= length(templateChoices)
            validity = "unidentified";
            choice = "";
            error = strcat("Invalid number of detected circles, should be ", num2str(length(templateChoices)), ", but was: ", num2str(length(ballotCircles)));
            return
        end
        
        %% - STEP 4
        %  - Figure out which circle(s) are marked
        markedCircleIndices = CheckMark(ballotCircles);
        
        %% - STEP 5
        %  - Figure out which circle is marked
        %  - If exactly one marked circle is found, declare the ballot valid and return 
        %  - otherwise, declare the ballot invalid and return
        if length(markedCircleIndices) == 1
            validity = "valid";
            choice = templateChoices(markedCircleIndices(1));
            return
        else
            validity = "invalid";
            choice = "";
            error = strcat("Invalid number of marked circles, should be 1, but was: ", num2str(length(markedCircleIndices)));
            % return all marked choices (if there is more than one)
            for i = markedCircleIndices
                choice = strcat(choice, templateChoices(i), " ");
            end
            return
        end
    catch e
        warning(getReport(e));
        validity = "invalid";
        error = "unknown (program error)";
        choice = "";
        return
    end
end

