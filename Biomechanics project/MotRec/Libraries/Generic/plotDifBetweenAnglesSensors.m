function plotDifBetweenAnglesSensors(FileName)

AngleMat = load([FileName,'_Angles']);
Angles = AngleMat.Angles;

SensorMat = load([FileName,'_Sensors']);
Sensors = SensorMat.Sensors;

Dif.angx1 = abs(Sensors.angx1 - Angles.angx1);
Dif.angy1 = abs(Sensors.angy1 - Angles.angy1);
Dif.angz1 = abs(Sensors.angz1 - Angles.angz1);
Dif.angx2 = abs(Sensors.angx2 - Angles.angx2);
Dif.angy2 = abs(Sensors.angy2 - Angles.angy2);
Dif.angz2 = abs(Sensors.angz2 - Angles.angz2);

Deltat = 0.01;
NFrames = size(Dif.angx1,1);
tend = Deltat * (NFrames-1);
time = 0.0:Deltat:tend;
Lim6 = 0.000001*ones(NFrames,1);
Lim12 = 0.000000000001*ones(NFrames,1);

figure('Name','Dif between Angle and Sensor in the axis X of Joint 1','NumberTitle','on','WindowStyle','docked'), axes('FontSize',20), hold on, grid on
% plot(time,Lim12,'r')
plot(time,Dif.angx1,'b')
% legend('Solver accuracy','Dif Angle Sensor',2) 
title({'Dif between Angle and Sensor in the axis X of Joint 1.';' Solver accuracy:1e-6 '},'FontSize',20)

% figure('Name','Dif between Angle and Sensor in the axis Y of Joint 1','NumberTitle','on','WindowStyle','docked'), axes('FontSize',20), hold on, grid on
% % plot(time,Lim12,'r')
% plot(time,Dif.angy1,'r')
% title(' axis Y of Joint 1. ','FontSize',20)
% 
% figure('Name','Dif between Angle and Sensor in the axis Z of Joint 1','NumberTitle','on','WindowStyle','docked'), axes('FontSize',20), hold on, grid on
% % plot(time,Lim12,'r')
% plot(time,Dif.angz1,'r')
% title(' axis Z of Joint 1.  ','FontSize',20)
% 
% figure('Name','Dif between Angle and Sensor in the axis X of Joint 2','NumberTitle','on','WindowStyle','docked'), axes('FontSize',20), hold on, grid on
% % plot(time,Lim12,'r')
% plot(time,Dif.angx2,'r')
% title(' axis X of Joint 2. ','FontSize',20)
% 
% figure('Name','Dif between Angle and Sensor in the axis Y of Joint 2','NumberTitle','on','WindowStyle','docked'), axes('FontSize',20), hold on, grid on
% % plot(time,Lim12,'r')
% plot(time,Dif.angy2,'r')
% title(' axis Y of Joint 2.','FontSize',20)
% 
% figure('Name','Dif between Angle and Sensor in the axis Z of Joint 2','NumberTitle','on','WindowStyle','docked'), axes('FontSize',20), hold on, grid on
% % plot(time,Lim12,'r')
% plot(time,Dif.angz2,'r')
% title('the axis Z of Joint 2. ','FontSize',20)