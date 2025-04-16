function  Phi = mkCtrMarkersNotInRef(BasisType, Points, Markers, BodyRef,SegName,InvR)

% MKCTRNOTINFEF creates constraints for the points not included in the reference
% as a function of the vectors of the reference
%
%   Phi = mkCtrPointsNotInRef(Basis, Points)
%   Inputs:
%     + BasisType is the type of Basis 
%     + Points is a vector of LocalPoints class
%     + Markers is a vector of LocalPoints class
%     + BodyRef is a struct (1x1) that contain:
%         Diri: strruc (1x1) that contain:
%               Global : sym value (3x1)
%               Local  : sym or double value (3x1)
%   Outputs:
%     + Phi are the symbolic constraints for points not included in the reference

% sizes
NMarkers = size(Markers,1);
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
for i=1:NMarkers
    MFin    = mkSymbolicXYZ(Markers(i).Point.Name);
    if isempty(Markers(i).LocCoord)
        MLocFin = mkSymbolicXYZ([SegName,'_',Markers(i).Point.Name]);
    else
        MLocFin = Markers(i).LocCoord;
    end
    LocU = MLocFin - PLocIni;
    U    = MFin - PIni;
    Phi = [Phi; mkCtrLinComb(U, LocU, BodyRef,InvR)];
end