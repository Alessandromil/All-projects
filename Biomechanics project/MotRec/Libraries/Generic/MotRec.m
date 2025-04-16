function MotRec(Path, FileXML)

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

% Define parameters
% Two options: define them in XML file or in M-File
% Here only XML file because this is for DHErgo Demonstrator interface
[Model,Experiment,Settings] = parseMotRecXML(Path, FileXML);

% Do reconstruction

MotionReconstruction(Model,Experiment,Settings);
