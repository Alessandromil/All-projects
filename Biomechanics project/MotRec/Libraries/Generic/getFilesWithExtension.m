function Files = getFilesWithExtension(Path, Ext, MessageON)

% getFileWithExtension
%    FileName = getFilesWithExtension(Path, Ext)
%
%  PATH: string containing path
%  EXT:  string containing file extension (without dot)
%  MessageON: 1 or 0, connects/disconnects messages
%
% This function searchs for files with extension EXT in folder PATH
%   - If there is only one file with extension EXT returns the filename
%   - If there is more than one file it prompts an Open File Dialog window
%     in order to choose one file with extension EXT
% It returns a file name with extension EXT

ListAll = dir(Path);
ListAll = ListAll(3:end);
Files   = {ListAll(~[ListAll.isdir]).name};
nFiles = length(Files);

FileNameExt = cell(nFiles,2);
for i=1:nFiles
    [FileNameExt{i,1}, FileNameExt{i,2}] = getFilenameAndExt(Files{i});   
end

IndexFiles = strmatch(Ext, FileNameExt(1:nFiles,2),'exact');
nFiles     = length(IndexFiles);

Files = cell(nFiles,1);

if nFiles == 0
    if MessageON
        disp(['Files with extension *.',Ext,' not found in folder ',Path]);
    end
else
    for i=1:nFiles
        Files{i} = [FileNameExt{IndexFiles(i),1}, '.', FileNameExt{IndexFiles(i),2}];
    end
end

