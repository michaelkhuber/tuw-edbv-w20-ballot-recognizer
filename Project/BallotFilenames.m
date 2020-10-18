function ballotFilenames = BallotFilenames()
    ballotFolder = "resources/ballots";
    
    if ~isfolder(ballotFolder)
        errorMessage = sprintf('Error: Ballot Folder does not exist:\n%s', ballotFolder);
        uiwait(warndlg(errorMessage));
        return;
    end
    
    filePattern = fullfile(ballotFolder, '*.png');
    ballotFiles = dir(filePattern);
    
    i = 1;
    ballotFilenames = strings(1, length(ballotFiles));
    for ballotFile = ballotFiles
        baseFilename = ballotFile.name;
        fullFilename = fullfile(ballotFolder, baseFilename);
        ballotFilenames(i) = fullFilename;
    end
end