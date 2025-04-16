function    Phi = mkCtrVectorsNotInRef(BasisType,Vectors,BodyRef,SegName,InvR)

% MKCTRVECTORSNOTINREF creates constraints for the vectors not included in the reference
% as a function of the vectors of the reference
%
%   Phi = mkCtrVectorsNotInRef(Vectors,SymAxes)
%   Inputs:
%   + BasisType is the type of Basis 
%   + Vectors is the vector of vector in segment
%   + SymAxes is a cell (3x2) that contain:
%       SymAxes{i,1}: Global symbolic coord.
%       SymAxes{i,2}: Local coord symbolic or not.

NVectors = size(Vectors,1);
% Check type of basis
% Basis formed by 3 vectors
if BasisType == 0
    IniVector = 4;
else
    IniVector = 3;
end
Phi = [];
for i=IniVector:NVectors
    V    = mkSymbolicXYZ(Vectors(i).Vector.Name);
    if isempty(Vectors(i).LocCoord)
        VLoc = mkSymbolicXYZ([SegName,'_',Vectors(i).Vector.Name]);
    else
        VLoc = Vectors(i).LocCoord;
    end
    Phi = [Phi; mkCtrLinComb(V, VLoc, BodyRef,InvR)];
end