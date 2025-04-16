function loadPath(varargin)

% LOADPATH adds to the matlab path the input paths
%   Inputs:
%     + Path1 is the first path that will be added to the matlab path - string
%     + Path2 is the second path that will be added to the matlab path - string
%       ...
%     + PathN is the N-th path that will be added to the matlab path - string
%   Outputs:
%     NONE

for i = 1: nargin
    Path = varargin{i};
    addpath(Path, '-begin');    
end
