%% ##########
%  ## MAIN ##
%  ##########
function ballotTable = Main()
% MAIN is the main function of our program. It reads *.jpg ballot images
% from the subfolder "resources/ballots" and tries to correctly find the voting
% choice that was marked for each one of them.
%
% Author:
%   Richard Binder
%
% Output:
%   ballotTable:    a table containing a row for each ballot image, with the
%   following columns:
%           index:          the index of each ballot image.
%           filename:     the file name of each ballot image.
%           success:     whether the marked choice(s) of the ballot image were correctly found 
%                               by our program or not. This is decided by comparing the found choices 
%                               to the expected choices given in the .csv file "resources/test_data.csv".
%           validity:        whether the vote on the ballot image is valid
%                               or invalid. (more or less than one marked choice makes the vote invalid). 
%                               This will be "unidentified" if the ballot image could not be processed correctly.
%           choices:      the choices that our program found for each ballot image.
%           errors:         if the validity of a ballot image is "invalid" or "unidentified", then
%                               the reason behind it is given here.

    % Add subfolders to path
    addpath(genpath(pwd));
    
    % Create results folder
    mkdir resources/results/Step2_Transform1;
    mkdir resources/results/Step2_Transform2;
    mkdir resources/results/Step3_Circles;

    % - Read in ordered ballot template choices
    templateChoices = Templ();

    % - Read in all Ballot Filenames from the Ballot Folder
    ballotFilenames = BallotFilenames();
    ballotIndices = 1:size(ballotFilenames,1);
    
    % - Manually choose Files (meant for debugging)
    ballotIndices = 17:size(ballotFilenames,1);
    
    % Get expected choices for each ballot from test_data.csv table
    testData = readtable('resources/test_data.csv');
    expectedChoices = string(testData.choices);
    
    % Only choose manually set ballots
    ballotFilenames = ballotFilenames(ballotIndices, :);
    expectedChoices = expectedChoices(ballotIndices);
    
    numBallots = size(ballotFilenames,1);
    % - Preallocate success array
    % - Each entry contains if the found choices are in line with the
    % expected choices from the test_data.csv table
    success = strings([numBallots, 1]);
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
        [success(i), validity(i), choices(i), errors(i)] = Pipeline(templateChoices, expectedChoices(i), ballotFilenames(i,:));
    end
    
    index = (1:numBallots)';
    colnames = ["index", "filename", "success", "validity", "choices", "errors"];
    names = ballotFilenames(:,2);
    ballotTable = table(index, names, success, validity, choices, errors, 'VariableNames', colnames);
    writetable(ballotTable, "resources/result.csv")
end

%% ##############
%  ## PIPELINE ##
%  ##############
function [success, validity, choice, error] = Pipeline(templateChoices, expectedChoice, ballotFilename)
%PIPELINE finds the marked choice(s) of a ballot image
%
% Author:
%   Richard Binder
%
% Input:
%   templateChoices:        the possible voting choices in the template
%   ballot.
%   expectedChoice:         the expected choice of the ballot image
%   ballotFilename:            the full path+name of the ballot image file
%
% Output:
%   success:        whether the marked choices were correctly found or not,
%   based on expectedChoice
%   validity:        whether the vote on the ballot image is valid
%   or invalid. (more or less than one marked choice makes the vote invalid). 
%   This will be "unidentified" if the ballot image could not be processed correctly.
%   choices:      the marked choices of the ballot image
%   errors:         if the validity of a ballot image is "invalid" or "unidentified", then
%   the reason behind it is given here.

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
        transformedBallot = Transform(ballotImg, ballotFilename(2), 1);
        %  - Transform2: Determine the table in the ballot paper, transform the image such
        %  that the ballot table is a (near) perfect rectangle. Then crop
        %  it to only the table.
        transformedBallot = Transform(transformedBallot, ballotFilename(2), 2);
        
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
        
        if strcmp(expectedChoice, choice)
            success = "true";
        else
            success = sprintf("false (Choice should be %s, but was %s)", expectedChoice, choice);
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

