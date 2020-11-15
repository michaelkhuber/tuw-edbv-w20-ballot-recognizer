function ballotImg = Read(ballotFilename)
    fprintf(1, 'Now reading ballot %s\n', ballotFilename);
    ballotImg = imread(ballotFilename);
end