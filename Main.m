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
    ballotFilenames = BallotFilenames()';

    % - Preallocate validity array
    % - Each entry can have take one of the three values:
    % - valid
    % - invalid
    % - unidentified
    validity = strings(size(ballotFilenames));
    % - Preallocate choices array
    % - Each entry contains the marked choices in a ballot (if any)
    choices = strings(size(ballotFilenames));
    
    % - Do something for every ballot image
    for i = 1:length(ballotFilenames)
        % Figure out the marked choice for a ballot image (if it is valid)
        [validity(i), choices(i)] = Pipeline(templateChoices, ballotFilenames(i));
    end
    
    index = (1:length(validity))';
    ballotTable = table(index, ballotFilenames, validity, choices);
    writetable(ballotTable, "result.csv")
end

%% ##############
%  ## PIPELINE ##
%  ##############
function [validity, choice] = Pipeline(templateChoices, ballotFilename)
        % - Implemented as suggested in the file
        % "Konzept_Wahlzettel_Erkennung.pdf" -> Point 5: Methodik
        % - For each Ballot, go through the following steps:
        
        %% - STEP 1
        %  - Read in ballot image
        ballotImg = Read(ballotFilename);
        %  - Prepare the image for circle matching
        preparedBallot = Prepare(ballotImg);
        
        %% - STEP 2
        %  - Match circles in ballot
        %  - output all circles in correct order
        ballotCircles = Circles(preparedBallot);
        
        %% - STEP 3
        %  - If the right amount of circles cannot be found, perform
        %  some transformation and try to match the circles again.
        if length(ballotCircles) ~= length(templateChoices)
            transformedBallot = Transform(preparedBallot, template);
            ballotCircles = Circles(transformedBallot);
        end
        %  - If the right amount of circles still cannot be found, declare the ballot's
        %  validity as unidentifiable -> cancel the pipeline and return
        if length(ballotCircles) ~= length(templateChoices)
            validity = "unidentified";
            choice = "";
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
            % return all marked choices (if there is more than one)
            for i = markedCircleIndices
                choice = strcat(choice, templateChoices(i), " ");
            end
            return
        end
end

