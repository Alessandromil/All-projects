function InverseDynamic(ExperimentPath,Model,Subject,Trial)
a=tic;

global PathBar

Index_barr= strfind(ExperimentPath,PathBar);
Experiment = ExperimentPath(Index_barr(end-1)+1:end-1);
RecPath = [ExperimentPath,'Subjects',PathBar,Subject,PathBar,'Results',PathBar];
RecName = [Model,'_',Experiment,'_',Trial,'_IK.mat'];
AceptedExtension{1,1} = 'mat';
checkFileAndPath(RecPath,RecName,AceptedExtension);
% RecName = 'PAM_MuscleLeftLeg_2Joint_Clutch_PedalINR1001_01_VS_CPR2_01_IK';
% ForcePath = [ExperimentPath,'\Subjects\',Subject,'\Motions\'];
% ForceFile = [Trial,'.for'];
RecMat = load([RecPath,RecName]);
Rec = RecMat.R;
Rec.ExperLogFileId = fopen([Rec.FileLogName(1:end-4),'_ID.log'],'w');
% ForceFileExist = exist([ForcePath,ForceFile],'file');
% if ForceFileExist ~=2
%     warning(['There is not force file for the trial', Trial])
% else
%     Rec.readForces(ForcePath,ForceFile);
% end
str = sprintf(['\n-------------------------------------------------------------------------\n', ...
                ' Inverse dynamic of model ',Model,' with trial ',Trial,...
                '\n-------------------------------------------------------------------------',]);
disp(str);
fprintf(Rec.ExperLogFileId,'%s', str);
Rec.doID(RecPath); 
[H, MI, S] = second2HMS(toc(a));
str =   ['Time: ',num2str(H),' hour(s) ',num2str(MI),' minute(s) and ',num2str(S),' second(s).'];
fprintf(Rec.ExperLogFileId,'%s\n', str);
fclose(Rec.ExperLogFileId);
disp(str);