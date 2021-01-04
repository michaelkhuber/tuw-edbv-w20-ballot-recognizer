function ballotImg = Read(ballotFilename)
% READ reads an image file based on its full file path+name
%
% Author:
%   Richard Binder
%
% Output:
%   ballotImg:      the image

    fprintf(1, 'Now reading ballot %s\n', ballotFilename);
    ballotImg = imread(ballotFilename);
end