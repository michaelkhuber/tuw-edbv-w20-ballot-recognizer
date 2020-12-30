function ballotFilenames = BallotFilenames()
    ballotFolder = "resources/ballots";
    
    % check if folder exists
    if ~isfolder(ballotFolder)
        errorMessage = sprintf('Error: Ballot Folder does not exist:\n%s', ballotFolder);
        uiwait(warndlg(errorMessage));
        return;
    end

    % Get existing files in folder
    filePattern = fullfile(ballotFolder, '*.jpg');
    ballotFiles = dir(filePattern);
    
    % Get Filenames from test_data.csv table
    testData = readtable('resources/test_data.csv');
    testDataNames = string(testData.filename)';
    
    % consistency check between test_data filenames and existing files
    for i = 1:length(ballotFiles)
        if( ~ismember(string(ballotFiles(i).name), testDataNames) )
            error("Inconsistent names: file %s was not found in test_data.csv table", string(ballotFiles(i).name));
        end
    end
    
    % Adding folder path to all filenames
    ballotFilenames = strings(length(testDataNames), 2);
    for i = 1:length(testDataNames)        
        ballotFilenames(i,1) = fullfile(ballotFolder, testDataNames(i));
        ballotFilenames(i,2) = testDataNames(i);
    end
end