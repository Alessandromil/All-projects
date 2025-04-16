function [IsValid, MessageError] = isValidPath(Path)

global PathBar PathBar2

PathInit = Path(1:3);

if ~strcmpi(Path(end),PathBar) % check is Path (full or relative) ends with the right character
    IsValid = 0; % invalid path
    MessageError = ['path must finish with "',PathBar2,'".'];
    
elseif strcmp(PathInit(1:2),PathBar2) || strcmp(PathInit(2:3),[':',PathBar]) || strcmp(PathInit(1:2),['.',PathBar])
    % IN WINDOWS
    %   '\\' -> Full path (network path),
    %   ':\') % Drive path (e.g. C:\myfolder\) Full path
    %   '.\' -> Path relative to xml file
    IsValid = 1; % valid path
    MessageError = [];
    
else
    MessageError = 'not valid path format';
    IsValid = 0; % invalid path
    
end
