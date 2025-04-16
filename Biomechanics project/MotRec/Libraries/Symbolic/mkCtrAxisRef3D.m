function Phi = mkCtrAxisRef3D(BodyRef)
%	MKCTR_AXISREF3D creates constraints for a 3D reference
% 
%   Phi = mkCtrAxisRef3D(SymAxes)
%   Inputs:    
%     + BodyRef is a struct (1x1) that contain:
%         Diri: strruc (1x1) that contain:
%               Global : sym value (3x1)
%               Local  : sym or double value (3x1)
%   Outputs:
%   + Phi are the symbolic constraints for the reference

Phi = [];
Phi = [Phi; mkCtrFixMod(BodyRef.Dir1.Global, BodyRef.Dir1.Local)];
Phi = [Phi; mkCtrFixMod(BodyRef.Dir2.Global, BodyRef.Dir2.Local)];
Phi = [Phi; mkCtrFixMod(BodyRef.Dir3.Global, BodyRef.Dir3.Local)];
Phi = [Phi; mkCtrFixAng(BodyRef.Dir1.Global, BodyRef.Dir2.Global, BodyRef.Dir1.Local, BodyRef.Dir2.Local)];
Phi = [Phi; mkCtrFixAng(BodyRef.Dir1.Global, BodyRef.Dir3.Global, BodyRef.Dir1.Local, BodyRef.Dir3.Local)];
Phi = [Phi; mkCtrFixAng(BodyRef.Dir2.Global, BodyRef.Dir3.Global, BodyRef.Dir2.Local, BodyRef.Dir3.Local)];

end

