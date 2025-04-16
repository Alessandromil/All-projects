%Indicate here all the motion files (.c3d) you wish to convert to .mat

% function C3D2MAT(InputPath,InputFile)
% 
% % receives a C3D file, extracts all data and saves part of it into a MAT file.
% % It is usefull when C3DServer software is not available. For example for
% % Apple computers

C3D2MAT('.\Experiments\ExpRightLeg\Motions\','G06_Ankle_DP.c3d')
C3D2MAT('.\Experiments\ExpRightLeg\Motions\','G06_Ankle_IE.c3d')
C3D2MAT('.\Experiments\ExpRightLeg\Motions\','G06_Knee_FE.c3d')
C3D2MAT('.\Experiments\ExpRightLeg\Motions\','G06_StaticPosture.c3d')



