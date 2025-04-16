function [VecNoNaN, indNaN] = getVecNoNaN(Vec)
% GETVECNONAN returns a struct with values and indices of vector Vec,
% which are not NaN. This is used to smooth a trajectory with gaps.
% The gaps are ignored and the trajectory is smoothed "piece by piece"
%
%  [VecNoNaN, indNaN] = getVecNoNaN(Vec)
%
%  Inputs:
%     Vec      : Array(nFrames x 1) that containts the trajectory of
%                only one coordinate.%
%                e.g. Vec=[NaN NaN 5 6 7 6 5 NaN NaN 3 4 5 3 NaN 2 2 1 NaN];
%  Outputs:
%      VecNoNaN: Struct with fields Val (values non NaN) & Id (indices)
%                For piece i of the trajectory we have
%                  VecNoNaN(i).Val
%                  VecNoNaN(i).Id
%      indNaN  : Struct with field Id (indices of NaNs)
%                For gap i of the trajectory we have
%                  indNaN(i).Id
%

% sizes
nFrames = length(Vec);
% initialize
cPieces = 0; % counter of non NaN values for current piece
nPieces = 1; % index of number of trajectory pieces
nGaps   = 0; % index of gaps

% logic vector with 1's in positions with NaN
VecNan = isnan(Vec);

% Check the first value of Vec
if VecNan(1) == 0 % if first value of Vec is not NaN
    cPieces = 1; % counter of non Nan values
    VecNoNaN(1).Val(cPieces) = Vec(1);
    VecNoNaN(1).Id(cPieces)  = 1;
else
    nGaps = 1;
    cNaNs = 1; % counter of NaN values for current gap.
    indNaN(nGaps).Id(cNaNs) = 1;
end

% Process the rest of the vector
for i=2:nFrames
    % IF previous value was non Nan and current value is NaN
    if VecNan(i-1) == 0 & VecNan(i) == 1
        nPieces = nPieces + 1; % increase counter of trajectory-piece
        nGaps   = nGaps  + 1; % increase counter of gaps
        cPieces = 0; % reset counter of non NaN values
        cNaNs   = 0; % reset counter of NaN values
    end
    
    % IF current value is non NaN
    if VecNan(i) == 0
        cPieces = cPieces+1;
        VecNoNaN(nPieces).Val(cPieces) = Vec(i);
        VecNoNaN(nPieces).Id(cPieces)  = i;
    else % IF current value is NaN
        cNaNs = cNaNs + 1;
        indNaN(nGaps).Id(cNaNs) = i;
    end
end

            
              
       
        
        
        
        
        
        
        
        
        
        