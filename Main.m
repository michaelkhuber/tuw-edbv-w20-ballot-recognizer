% - Our pipeline can return 4 different outputs:
% - ballot is invalid
% - ballot's validity is unidentifiable
% - ballot is inaproppriately marked
% - ballot is valid and apropproiately marked -> return marked choice


%% ##########
%  ## MAIN ##
%  ##########
function ballotTable = Main()
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
    ballotIndices = 1:94;
    
    % Get expected choices for each ballot from test_data.csv table
    testData = readtable('resources/test_data.csv');
    testDataChoices = string(testData.choices);
    
    % Only choose manually set ballots
    ballotFilenames = ballotFilenames(ballotIndices, :);
    testDataChoices = testDataChoices(ballotIndices);
    
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
        [success(i), validity(i), choices(i), errors(i)] = Pipeline(templateChoices, testDataChoices(i), ballotFilenames(i,:));
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

