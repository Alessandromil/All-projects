function calcSubPar(ExperimentPath, SubParFiles)

global PathBar

NSubParFiles = size(SubParFiles,2);
for i=1:NSubParFiles
    File = SubParFiles{i};
    FilePath = File(1:(end-4));
    Path = [ExperimentPath,'Subjects',PathBar,FilePath,PathBar];
    AceptedExtension{1,1} = 'alm';
    checkFileAndPath(Path,File,AceptedExtension);    
    RealHuman = HUMAN_MODEL();
    RealHuman.parseAlmFile(Path,File);
    RealHuman.addModelJoints();
    RealHuman.calcTCS();
    RealHuman.calcLCS();
    RealHuman.calcInertiaPars();
    RealHuman.writeSubjectParsIntermed(File(1:(end-4)),Path);
%     RealHuman.drawLCS();
%   ramsispar
    RealHuman.calcALLCS_T_RamsisLCS();
    RealHuman.writeRSP(File(1:(end-4)),Path);
    % RealHuman.pampar(ParamentrosPam);
end
