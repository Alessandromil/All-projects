function C3D2MAT(InputPath,InputFile)

% receives a C3D file, extracts all data and saves part of it into a MAT file.
% It is usefull when C3DServer software is not available. For example for
% Apple computers


% read C3D file
[Markers, Frequency, IndexAllVisible, MarkerNames, MarkerCoords] = readC3D(InputPath, InputFile);

% get only filename
[Filename, ~] = getFilenameAndExt(InputFile);                

% save data in a MAT file.
disp(['Saving motion data in file ',InputFile,' ...'])
save([InputPath,Filename],'Frequency','MarkerNames','MarkerCoords');
disp(['Process finished!'])
