function FilesWithPattern = getFilesWithPattern(Path, Prefix, Suffix)

% getFilesWithPatern
%    Files = getFilesWithPatern(Path, Prefix, Sufix, MessageON)
%
%  PATH  : string containing path
%  PREFIX: string containing file prefix 
%  SUFFIX: string containing file suffix 
%  MessageON: 1 or 0, connects/disconnects messages
%
% This function searchs for files with pattern [PREFIX,*,SUFFIX]
% It returns a cell with all files with pattern

ListAll = dir(Path);
ListAll = ListAll(3:end);
FilesAll = {ListAll(~[ListAll.isdir]).name};

% Get files starting by string Prefix
IndexFiles = strmatch(Prefix,FilesAll);
FilesWithPrefix = FilesAll(IndexFiles)';


% Flip files in order to check suffix with function strmatch
nFilesWithPrefix = length(FilesWithPrefix);
FilesWithPrefix_fliplr = cell(nFilesWithPrefix,1);
for i=1:length(FilesWithPrefix)
    FilesWithPrefix_fliplr{i} = fliplr(FilesWithPrefix{i});
end
Suffix_fliplr = fliplr(Suffix);


% Get files ending by string Suffix
IndexFiles = strmatch(Suffix_fliplr,FilesWithPrefix_fliplr);
FilesWithPattern = FilesWithPrefix(IndexFiles);


