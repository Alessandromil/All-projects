function SegRef = mkSymbSegRef(BasisType,Points,Vectors,SegName)

% MKSYMSEGREF estimated symbolic axes depending of type of basis
%
%   SegRef = mkSymbSegRef(BasisType,Points,Vectors)
%   Inputs:
%   + Basis is a class and contain type of basis
%   + Points is a vector of LocalPoint class in segments
%   + Vector is a vector of LocalVector class in segments
%   Outputs:
%   + SegRef is a struct (1x1) that contain:
%       Diri: strruc (1x1) that contain:
%             Global : sym value (3x1)
%             Local  : sym or double value (3x1)

% Put in symbolic form
Axis1 = mkSymbolicXYZ(Vectors(1).Vector.Name);
Axis2 = mkSymbolicXYZ(Vectors(2).Vector.Name);
if isempty(Vectors(1).LocCoord)
    LocAxis1 = mkSymbolicXYZ([SegName,'_',Vectors(1).Vector.Name]);
else
    LocAxis1 = Vectors(1).LocCoord;
end
if isempty(Vectors(2).LocCoord)
    LocAxis2 = mkSymbolicXYZ([SegName,'_',Vectors(2).Vector.Name]);
else
    LocAxis2 = Vectors(2).LocCoord;
end

% Check type of basis
% Basis formed by 3 vectors
if BasisType == 0
    Axis3 = mkSymbolicXYZ(Vectors(3).Vector.Name);
    if isempty(Vectors(3).LocCoord)
        LocAxis3 = mkSymbolicXYZ([SegName,'_',Vectors(3).Vector.Name]);
    else
        LocAxis3 = Vectors(3).LocCoord;
    end
% Basis formed by 2 vectors and 2 points    
elseif BasisType == 1
    P1 = mkSymbolicXYZ(Points(1).Point.Name);
    P2 = mkSymbolicXYZ(Points(2).Point.Name);
    Axis3 = P2-P1;
    if isempty(Points(1).LocCoord)
        LocP1 = mkSymbolicXYZ([SegName,'_',Points(1).Point.Name]);
    else
        LocP1 = Points(1).LocCoord;
    end
    if isempty(Points(2).LocCoord)
        LocP2 = mkSymbolicXYZ([SegName,'_',Points(2).Point.Name]);
    else
        LocP2 = Points(2).LocCoord;
    end
    LocAxis3 = LocP2-LocP1;
else
    error('Incorrect basis')
end

SegRef.Dir1.Global = Axis1;
SegRef.Dir1.Local  = LocAxis1;
SegRef.Dir2.Global = Axis2;
SegRef.Dir2.Local  = LocAxis2;
SegRef.Dir3.Global = Axis3;
SegRef.Dir3.Local  = LocAxis3;
