function varargout = checkAndCorrectFileAndPath(Path,File,AceptedExtension)
% CHECKFILEANDPATH checks if the extension is correct and if exist the path and the file. 
%   checkFileAndPath(Path,File,AceptedExtension)
%   Inputs:
%       + Path                            string 
%       + File                            string
%       + AceptedExtension                cell{NExtensionsx1}
%         There are all the posible extensions for this type of file

global PathBar

% 1) Check the acepted extensions
[FileName , Extension] = getFilenameAndExt(File);
NExtension = size(AceptedExtension,1);
ContExten = 0;
for i=1:NExtension
    if strcmpi(Extension,AceptedExtension{i})
        ContExten = 1;
        break;
    end
end
if ContExten == 0
    AllExtensions = [];
    for i=1:NExtension
        AllExtensions = [AllExtensions,'.',AceptedExtension{i},' or '];
    end
    AllExtensions(end-3:end)= '.   ';
    if isempty(FileName)
        error(['File "',Extension,'" -> acepted extensions are only:',AllExtensions]);
    else
        error(['File "',FileName,'.',Extension,'" -> acepted extensions are only:',AllExtensions]);
    end
end
% 2) Check the path
if ~strcmpi(Path(end),PathBar)
    Path = [Path,PathBar];    
end
PathExist = exist(Path,'dir');
if PathExist~=7
    error(['The path: "',Path,'" does not exist.'])
end
% 3) Check the file
FileExist = exist([Path,File],'file');
if FileExist ~=2
    error(['File ',File,' can not be found in path ',getPrintPath(Path)])
end


% define output vars if necessary
if nargout == 1
    varargout{1} = Path;
end