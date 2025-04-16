function loadBioToolboxPath(BioToolboxPath)

% LOADBIOTOOLBOXPATH adds to the Matlab path the path of the 
% Biomechanics Toolbox which are:
%   - Parameters Library Path 

%   Inputs:
%     + BioToolboxPath is the path where biomechanis Toolbox is installed - string
%       E.g. 'D:\Users\BioToolbox' (Windows)
%   Outputs:
%     NONE

global PathBar

% ------------------------------------------------------------------------------------------------
% Add the Libraries Path 
% ------------------------------------------------------------------------------------------------
% Make model Library
addpath([BioToolboxPath,PathBar,'Libraries',PathBar,'Classes'], '-begin');
% Make symbolic Library
addpath([BioToolboxPath,PathBar,'Libraries',PathBar,'Symbolic'], '-begin');
% Generic Library
addpath([BioToolboxPath,PathBar,'Libraries',PathBar,'Generic'],'-begin');
% C3D 
addpath([BioToolboxPath,PathBar,'Libraries',PathBar,'readC3D'],'-begin');

