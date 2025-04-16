function  Phi = mkCtrPointsNotInRef(BasisType, Points, BodyRef,SegName,InvR)

% MKCTRNOTINFEF creates constraints for the points not included in the reference
% as a function of the vectors of the reference
%
%   Phi = mkCtrPointsNotInRef(Basis, Points)
%   Inputs:
%     + BasisType is the type of Basis 
%     + Points is a vector of LocalVector class
%     + BodyRef is a struct (1x1) that contain:
%         Diri: strruc (1x1) that contain:
%               Global : sym value (3x1)
%               Local  : sym or double value (3x1)
%   Outputs:
%     + Phi are the symbolic constraints for points not included in the reference

% sizes
NPoints = size(Points,1);
% Check type of basis
% Basis formed by 3 vectors
if BasisType == 0
    IniPoint = 2;
else
    IniPoint = 3;
end
PIni = mkSymbolicXYZ(Points(1).Point.Name);
if isempty(Points(1).LocCoord)
        PLocIni = mkSymbolicXYZ([SegName,'_',Points(1).Point.Name]);
else
        PLocIni = Points(1).LocCoord;
end
Phi =[];
for i=IniPoint:NPoints
    PFin    = mkSymbolicXYZ(Points(i).Point.Name);
    if isempty(Points(1).LocCoord)
        PLocFin = mkSymbolicXYZ([SegName,'_',Points(i).Point.Name]);
    else
        PLocFin = Points(i).LocCoord;
    end
    LocU = PLocFin - PLocIni;
    U    = PFin - PIni;
    Phi = [Phi; mkCtrLinComb(U, LocU, BodyRef,InvR)];
end
