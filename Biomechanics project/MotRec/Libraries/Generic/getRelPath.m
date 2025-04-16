function RelPath = getRelPath(ReferencePath, OtherPath)

% RELPATH returns the relative path to go from ReferencePath to OtherPath
%
%   RelPath = getRelPath(ReferencePath, OtherPath)
%   Inputs:
%     + ReferencePath is the refence path - string
%     + OtherPath is the objective path - string
%   Outputs:
%     + RelPath is the relative path to go from the ReferencePath to the
%       OtherPath


% directories to test the function
% PlaybackPath = 'D:\Users\sausejo\proyectos\2134-ErgoPSA\Development\Experiments\braking\Exp1\resultsIK\';
% GraphicsPath = 'D:\Users\sausejo\proyectos\2134-ErgoPSA\Development\BioModels\InvKinModels\RLowerLimb\';
%
% PlaybackPath = '\p\';
% GraphicsPath = '\p\';
%
% PlaybackPath = '\p\a\b\';
% GraphicsPath = '\p\';
%
% PlaybackPath = '\p\';
% GraphicsPath = '\p\a\b\';

global PathBar

% find the length of the shortest string
LenghReferPath = length(ReferencePath);
LenghOtherPath = length(OtherPath);
if LenghReferPath >= LenghOtherPath
    imax = LenghOtherPath;
else
    imax = LenghReferPath;
end

% find the location where the two paths are diferent
i = 1;
while strncmp(ReferencePath(1:i), OtherPath(1:i), i)
    i = i + 1;
    if i > imax
        break
    end
end

% Store only the part of the paths that are different
OtherPath_RelPath     = OtherPath(i:end);
ReferencePath_RelPath = ReferencePath(i:end);

% count the number of directories that we have to go up from
% the ReferencePath until the OtherPath.
NumberDirs = 1;
[Dir, Reminder] = strtok(ReferencePath_RelPath,PathBar);
if isempty(Dir)
        % we do not need to go up
    DirsUp = ['.',PathBar];
    
else
    while ~strcmp(Reminder,PathBar) 
        [Dir, Reminder] = strtok(Reminder,PathBar);
        NumberDirs = NumberDirs + 1;
    end
    
    % create the string to go up NumberDirs directories
    DirsUp = '';
    for i = 1 : NumberDirs
        DirsUp = [DirsUp, '..',PathBar];
    end
    
end

% relative path between ReferencePath and OtherPath
RelPath = [DirsUp,OtherPath_RelPath];
