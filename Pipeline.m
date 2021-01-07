%% ##############
%  ## PIPELINE ##
%  ##############
function [success, validity, choice, error] = Pipeline(templateChoices, testDataChoice, ballotFilename)
    error = "";
    try
        % - Implemented as suggested in our concept file
        % "Konzept_Wahlzettel_Erkennung.pdf" -> Point 5: Methodik
        % - For each Ballot, go through the following steps:
        
        %% - STEP 1
        %  - Read in ballot image
        ballotImg = Read(ballotFilename(1));
        
        %% - STEP 2
        %  - Transform1: Determine the ballot paper in the image, transform the image such
        %  that the ballot paper is a (near) perfect rectangle. Then crop
        %  it to only the ballot paper.
        [transformedBallot, step] = Transform(ballotImg, ballotFilename(2), 1);
        %  - Transform2: Determine the table in the ballot paper, transform the image such
        %  that the ballot table is a (near) perfect rectangle. Then crop
        %  it to only the table.
        if step== 1
            transformedBallot = Transform(transformedBallot, ballotFilename(2), 2);
        end
        
        %% - STEP 3
        %  - Match circles in ballot
        %  - output all circles in correct order
        ballotCircles = Circles(transformedBallot, ballotFilename(2));
        
        %  - If the right amount of circles cannot be found, declare the ballot's
        %  validity as unidentified -> cancel the pipeline and return
        if length(ballotCircles) ~= length(templateChoices)
            success = "false (Invalid number of detected circles)";
            validity = "unidentified";
            choice = "";
            error = strcat("Invalid number of detected circles, should be ", num2str(length(templateChoices)), ", but was: ", num2str(length(ballotCircles)));
            return
        end
        
        %% - STEP 4
        %  - Figure out which circle(s) are marked
        markedCircleIndices = CheckMark(ballotCircles);
        
        %% - STEP 5
        %  - If exactly one marked circle is found, declare the ballot valid and return 
        %  - otherwise, declare the ballot invalid and return
        if length(markedCircleIndices) == 1
            validity = "valid";
            choice = templateChoices(markedCircleIndices(1));
        else
            validity = "invalid";
            choice = "";
            error = strcat("Invalid number of marked circles, should be 1, but was: ", num2str(length(markedCircleIndices)));
            % return all marked choices (if there is more than one)
            for i = markedCircleIndices
                choice = strcat(choice, templateChoices(i));
                if i ~= markedCircleIndices(length(markedCircleIndices))
                    choice = strcat(choice, " ");
                end
            end
        end
        
        if strcmp(testDataChoice, choice)
            success = "true";
        else
            success = sprintf("false (Choice should be %s, but was %s)", testDataChoice, choice);
        end
        return
    catch e
        warning(getReport(e));
        success = "false (program error)";
        validity = "unidentified";
        error = "unknown (program error)";
        choice = "";
        return
    end
end