function dispWarning(MessageLines, InfoType, varargin)

% MessageLines  string containing warning message
% InfoType      integer, 0 minimun info, 1 standard info
% varargin{1}   FileID, file identifier in which the warning must also be written   

SepTop    = '\n=========================================================================================\n';
SepBottom = '\n=========================================================================================\n';

if ~iscell(MessageLines)
    MessageLines = {MessageLines};
end
nLines = length(MessageLines);
Gap = '           ';

% For minimun information
Message = ['  WARNING: ',MessageLines{1}];
for i=2:nLines
    Message = [Message,Gap,MessageLines{i},'\n'];    
end
% For full information
FullMessage = sprintf([SepTop,Message,SepBottom]);

% disp data in command window
if InfoType == 0
    disp(sprintf(Message));
elseif InfoType == 1
    disp(FullMessage);
end

% write in file
if nargin == 3
    FileID = varargin{1};
    fprintf(FileID, '%s\n', FullMessage);
end
