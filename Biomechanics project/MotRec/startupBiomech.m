% STARTUPLIBRARY Startup file for the Biomechanics Toolbox
% The file is executed when MATLAB starts up, if it exists:
%   - Anywhere on the path
%   - 'Start in' directory 

% ------------------------------------------------------------------------
% MAC or Window?
% ------------------------------------------------------------------------

global PathBar PathBar2
if ismac
    % Code to run on Mac platform
    PathBar  = '/';
    PathBar2 = '/';
elseif ispc
    % Code to run on Windows platform
    PathBar  = '\';
    PathBar2 = '\\'; 
end

% ------------------------------------------------------------------------
% Toolbox settings
% ------------------------------------------------------------------------
% define Biomechanics Toolbox PATH
BioToolboxPATH = pwd;    % get actual directory

% add Biomechanics Toolbox Paths to Matlab path
% add Library classes to path
addpath([BioToolboxPATH,PathBar,'Libraries',PathBar,'Classes'], '-begin');
% add Generic library to path
addpath([BioToolboxPATH,PathBar,'Libraries',PathBar,'Generic'],'-begin');
% add Symbolic custom library to path
addpath([BioToolboxPATH,PathBar,'Libraries',PathBar,'readC3D'], '-begin');
% add Symbolic custom library to path
addpath([BioToolboxPATH,PathBar,'Libraries',PathBar,'Symbolic'], '-begin');


% display info
disp('Biomechanics Toolbox is installed:');
disp(['  Location   : ',BioToolboxPATH]);
disp('  Description: Creates and simulates biomechanics models.');

% ------------------------------------------------------------------------
% Work Path settings
% ------------------------------------------------------------------------
% define Work directory PATH
WorkPATH = BioToolboxPATH;

% Set Matlab's current directory to our Work directory
cd(WorkPATH);

% delete variables
clear BioToolboxPATH WorkPATH PathBar PathBar2
