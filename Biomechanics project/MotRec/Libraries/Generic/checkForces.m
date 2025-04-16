function checkForces(FilePath,FileName)

Sub = load([FilePath,FileName]);
mt = 0;
for i=2:size(Sub.SaveSubject.Segments,1)
    mi = Sub.SaveSubject.Segments(i).Mass;
    mt = mt+mi;
end
NFrames = size(Sub.SaveSubject.Segments(2).Proximal.Joint.F,1);
Weigh = [0;0;-9.81]*mt;
for i=1:NFrames
    NormFPelvis(i) = norm(Sub.SaveSubject.Segments(2).Proximal.Joint.F(i,1:3));
    FExt_We = Sub.SaveSubject.Segments(5).F_Ext.Value(i,1:3)'+ Weigh;
    NormFExt_We(i) = norm(FExt_We);
    Dif(i) = NormFPelvis(i)-NormFExt_We(i);
end
tend = 0.01 * (NFrames-1);
time = 0.0:0.01:tend;
Figurei = figure('Name',FileName(1:end-4),'NumberTitle','off','WindowStyle','docked');
axes('FontSize',12);
hold on
grid on
plot(time,NormFPelvis,'r');
plot(time,NormFExt_We,'b');
saveas(Figurei,FileName(1:end-4))