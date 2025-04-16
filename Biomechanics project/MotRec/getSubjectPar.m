
StaticPath = [pwd,'\Experiments\ExpRightLeg\Motions\']; 
StaticFile = 'Group10_static1.c3d'; 

[Markers, Frequency, IndexAV, MarkerNames, MarkerCoords]=readC3D(StaticPath, StaticFile);


Glob_G10_RT = Markers.MarkerSet_G10_RT(IndexAV,:)';
Glob_G10_LFE = Markers.MarkerSet_G10_LFE(IndexAV,:)';
Glob_G10_MFE = Markers.MarkerSet_G10_MFE(IndexAV,:)';
Glob_G10_RS = Markers.MarkerSet_G10_RS(IndexAV,:)';
Glob_G10_LM = Markers.MarkerSet_G10_LM(IndexAV,:)';
Glob_G10_MM = Markers.MarkerSet_G10_MM(IndexAV,:)';
Glob_G10_F1 = Markers.MarkerSet_G10_F1(IndexAV,:)';
Glob_G10_F2 = Markers.MarkerSet_G10_F2(IndexAV,:)';
Glob_G10_F3 = Markers.MarkerSet_G10_F3(IndexAV,:)';
 
Glob_AJC = (Glob_G10_LM+Glob_G10_MM)/2;
Glob_KJC = (Glob_G10_LFE+Glob_G10_MFE)/2;
Glob_HJC = Glob_KJC+[0;0.47;0];


% Shank
Glob_Os = Glob_AJC;
bs = Glob_KJC-Glob_AJC;
as = Glob_G10_LFE-Glob_G10_MFE;
Glob_Xs = cross(bs,as)/(norm(cross(bs,as)));
Glob_Ys = (Glob_KJC-Glob_AJC)/(norm(Glob_KJC-Glob_AJC));
Glob_Zs = cross(Glob_Xs,Glob_Ys);

% Foot
Glob_Of = Glob_AJC;
Glob_Yf = Glob_Ys;
Glob_Zf = Glob_Zs;
Glob_Xf = Glob_Xs;

%Thigh 
Glob_Ot = Glob_KJC;
Glob_Yt = (Glob_HJC-Glob_KJC)/(norm(Glob_HJC-Glob_KJC));
at=as/norm(as);
Glob_Xt = cross(Glob_Yt,at);
Glob_Zt = cross(Glob_Xt,Glob_Yt);



% Shank --------------------
Glob_R_Shank = [Glob_Xs,Glob_Ys,Glob_Zs];

Shank_RS = (Glob_R_Shank')*(Glob_G10_RS-Glob_Os);
Shank_LM = (Glob_R_Shank')*(Glob_G10_LM-Glob_Os);
Shank_MM = (Glob_R_Shank')*(Glob_G10_MM-Glob_Os);
Shank_AJC = (Glob_R_Shank')*(Glob_AJC-Glob_Os);
Shank_KJC = (Glob_R_Shank')*(Glob_KJC-Glob_Os);

ShankData.Name = 'Shank';
ShankData.MarkerNames  = {'RS', 'LM', 'MM', 'AJC', 'KJC'};
ShankData.MarkerCoords = [Shank_RS, Shank_LM, Shank_MM, Shank_AJC, Shank_KJC];      
drawLCS(ShankData,        1,          'Shank');



% Foot ---------------------
Glob_R_Foot = [Glob_Xf,Glob_Yf,Glob_Zf];

Foot_F1 = (Glob_R_Foot')*(Glob_G10_F1-Glob_Of);
Foot_F2 = (Glob_R_Foot')*(Glob_G10_F2-Glob_Of);
Foot_F3 = (Glob_R_Foot')*(Glob_G10_F3-Glob_Of);
Foot_AJC = (Glob_R_Foot')*(Glob_AJC-Glob_Of);
Foot_KJC = (Glob_R_Foot')*(Glob_KJC-Glob_Of);


FootData.Name = 'Foot';
FootData.MarkerNames  = {'F1', 'F2', 'F3', 'AJC'};
FootData.MarkerCoords = [Foot_F1, Foot_F2, Foot_F3, Foot_AJC];      
drawLCS(FootData,        1,          'Foot');



% Thigh --------------------
Glob_R_Thigh = [Glob_Xt,Glob_Yt,Glob_Zt];

Thigh_RT = (Glob_R_Thigh')*(Glob_G10_RT-Glob_Ot);
Thigh_LFE = (Glob_R_Thigh')*(Glob_G10_LFE-Glob_Ot);
Thigh_MFE = (Glob_R_Thigh')*(Glob_G10_MFE-Glob_Ot);
Thigh_HJC = (Glob_R_Thigh')*(Glob_HJC-Glob_Ot);
Thigh_KJC = (Glob_R_Thigh')*(Glob_KJC-Glob_Ot);

ThighData.Name = 'Thigh';
ThighData.MarkerNames  = {'RT', 'LFE', 'MFE', 'HJC', 'KJC'};
ThighData.MarkerCoords = [Thigh_RT, Thigh_LFE, Thigh_MFE, Thigh_HJC, Thigh_KJC];      
drawLCS(ThighData,        1,          'Thigh');




