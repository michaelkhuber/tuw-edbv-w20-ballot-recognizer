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
    ballotIndices = 1:110;
    
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

