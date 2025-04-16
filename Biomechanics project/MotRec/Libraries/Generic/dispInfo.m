function dispInfo(InfoType, FullMessage, varargin)
%function dispInfo(InfoType, FullMessage, MinMessage, FileID)


if nargin == 3
    if isnumeric(varargin{1})
        FileID = varargin{1};
        fprintf(FileID, '%s\n', FullMessage);
    else
        MinMessage = varargin{1};
    end
elseif nargin == 4
    MinMessage = varargin{1};
    FileID = varargin{2};
    fprintf(FileID, '%s\n', FullMessage);
end

if(InfoType == 1)
    disp(FullMessage);
elseif (InfoType == 0) && exist('MinMessage','var') == 1
    disp(MinMessage);
end