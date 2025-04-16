function [SmarkIndex, nSmarkInVecCoord] = getSmarkIndex(VecCoord, SmarkNames)

% GETSMARKINDEX returns a vector with the index of each skin-marker's
% X-coordinate location in q or z. Variable 'VecCoord' can be q or z
%
%   [SmarkIndex, nSmarkInVecCoord] = getSmarkIndex(VecCoord, SmarkNames)
%   Inputs:
%     + SmarkNames is a cell (nSmark x 1) with the names of the skin-markers.
%       Each element of the cell is a string.
%     + VecCoord is a symbolic array(nCoord x 1) with symbolic
%       variables. It can be q (vector of generalized coor.) or
%       z (vector of inputs).
%   Outputs:
%     + SmarkIndex is the integer array with the indexes of each
%       skin-marker's X-coordinate location in symbolic vector 'VecCoord'.
%     + nSmarkInVecCoord is the number of skin-markers in the symbolic
%       vector 'VecCoord'.


% sizes
nCoord      = length(VecCoord);
nSmarkNames = length(SmarkNames);

% initialize
InputIndex = 1;
nSmarkInVecCoord = 0;
SmarkIndex = [];

while InputIndex <= nCoord
    InputName = char(VecCoord(InputIndex));
    InputName = InputName(1:end-1);
        
    for i=1:nSmarkNames
        if strcmp(SmarkNames{i},InputName)
            nSmarkInVecCoord = nSmarkInVecCoord + 1;
            SmarkIndex(nSmarkInVecCoord) = InputIndex;
            InputIndex = InputIndex + 2;
            break
        end
    end
    InputIndex = InputIndex + 1;  
end
