function PrintPath = getPrintPath(Path)

% GETPRINTPATH transforms the string 'Path' that contains a path in
% a string that contains a printable version of Path. This is required
% because functions like FPRINTF or SPRINTF need a special format to be
% able to print a path appropriatedly. To print a backslash character 
% or percent character two blackslash or two percent characters are 
% needed: \\ or  %% 
%  
%   PrintPath = getPrintPath(Path) 
%   Inputs:
%     + Path is a string with a path, e.g. 'D:\User\Work'
%   Outputs:
%     + PrintPath is a string with a path that can be printed with
%       functions FPRINTF or SPRINTF, e.g. 'D:\\User\\Work'
global PathBar PathBar2

PrintPath = strrep(Path,PathBar,PathBar2);
