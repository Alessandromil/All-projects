function Experiment = addDataFiles(varargin)

% Extract input data from variable varargin
if nargin == 2 || nargin == 3    
    Experiment  = varargin{1};
    MotionFiles = varargin{2};
end
if nargin == 3    
    SubDataFiles = varargin{3};
end

% do things
if nargin == 2
    TmpDataFiles{1,1} = [];
    TmpDataFiles{1,2} = MotionFiles;
elseif nargin == 3
    TmpDataFiles{1,1} = SubDataFiles;
    TmpDataFiles{1,2} = MotionFiles;    
end
if ~isfield(Experiment,'DataFiles')
    Experiment.DataFiles = {};
end
Experiment.DataFiles =[Experiment.DataFiles;TmpDataFiles];