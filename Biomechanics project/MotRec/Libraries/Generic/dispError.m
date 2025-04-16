function dispError(Message, InfoType, varargin)


SepTop    = '\n=========================================================================================\n';
SepBottom = '\n=========================================================================================';

FullMessage = sprintf([SepTop,' ERROR: ',Message,SepBottom]);

if nargin == 3
    FileID = varargin{1};
    fprintf(FileID, '%s\n', FullMessage);
end
if InfoType == 0 || InfoType == 1
    error(FullMessage);
end

