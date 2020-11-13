function ballotFilenames = BallotFilenames()
    ballotFolder = "resources/ballots";
    
    if ~isfolder(ballotFolder)
        errorMessage = sprintf('Error: Ballot Folder does not exist:\n%s', ballotFolder);
        uiwait(warndlg(errorMessage));
        return;
    end
    
    filePattern = fullfile(ballotFolder, '*.png');
    ballotFiles = dir(filePattern);
    bla1 = ballotFiles(1);
    bla2 = ballotFiles(2);
    bla3 = ballotFiles(3);
    
    ballotFilenames = strings(1, length(ballotFiles));
    for i = 1:length(ballotFiles)
        ballotFilenames(i) = fullfile(ballotFolder, ballotFiles(i).name);
    end
end