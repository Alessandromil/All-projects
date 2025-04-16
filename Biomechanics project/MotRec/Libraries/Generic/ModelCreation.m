function ModelCreation(Models)
%ModelCreation Summary of this function goes here
%   Detailed explanation goes here
NModels = size(Models,1);

for i=1:NModels
    c=tic;
    ModelPath = Models(i).Path;
    ModelFile = Models(i).File;
    AceptedExtension{1,1} = 'xml';
    checkFileAndPath(ModelPath,ModelFile,AceptedExtension);
    disp(' ');
    disp(['Creating model ',ModelFile(1:end-4)]);
    disp(['  in folder ',ModelPath])
    Human = HUMAN_MODEL();
    Human.parseXmlModelFile(ModelPath,ModelFile);
    Human.fillq;
    Human.addMarkerInq();
    Human.mkModelCtrs();
    save([ModelPath,ModelFile(1:end-4)],'Human');
    [H, MI, S] = second2HMS(toc(c));
    printElapsedTime(H, MI, S, 3, 'Model creation time: ');
end

