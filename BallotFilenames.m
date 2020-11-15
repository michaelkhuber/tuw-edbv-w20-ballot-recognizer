function ballotFilenames = BallotFilenames()
    ballotFolder = "resources/ballots";
    
    if ~isfolder(ballotFolder)
        errorMessage = sprintf('Error: Ballot Folder does not exist:\n%s', ballotFolder);
        uiwait(warndlg(errorMessage));
        return;
    end
    
    testData = readtable('resources/test_data.csv');
    fileNames = testData.filename;

    filePattern = fullfile(ballotFolder, '*.jpg');
    ballotFiles = dir(filePattern);
    bla1 = ballotFiles(1);
    bla2 = ballotFiles(2);
    bla3 = ballotFiles(3);
    
    ballotFilenames = strings(1, length(fileNames));
    for i = 1:length(ballotFiles)
        ballotFilenames(i) = fullfile(ballotFolder, fileNames(i));
    end
end