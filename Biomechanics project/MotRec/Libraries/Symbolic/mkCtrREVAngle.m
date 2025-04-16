function Phi = mkCtrREVAngle(Angle)

% MKCTRREVANGLE creates the relative angle constraints for a Revolution joint
%
%   Phi = mkCtrREVAngle(Angle)
%   Inputs:

%   Outputs:
%     + Phi are the symbolic Revolute joint constraints
% set variables
Axe1  = Angle.Joint.Vector1.Vector.Name; % Axis of the Revolution joint. It MUST be a UNITARY Vector
VAxe1Seg1 = Angle.V1Seg1.Vector.Name; % Reference vector in Seg1. It MUST be a UNITARY Vector
VAxe1Seg2 = Angle.V1Seg2.Vector.Name; % Reference vector in Seg2. It MUST be a UNITARY Vector
Alpha = Angle.Name1; % Relative angle of the REVJoint


% Check each input. If it is a char variable create the corresponding symbolic variable
if ischar(Axe1), Glob_Vec_Axe1 = mkSymbolicXYZ(Axe1);  end
if ischar(VAxe1Seg1), Glob_Vec_VAxe1Seg1 = mkSymbolicXYZ(VAxe1Seg1);  end
if ischar(VAxe1Seg2), Glob_Vec_VAxe1Seg2 = mkSymbolicXYZ(VAxe1Seg2);  end
if ischar(Alpha), Alpha = mkSymbAngle(Alpha);  end

% -----------------------------------------------
% Rotation
% -----------------------------------------------

Glob_U   = Glob_Vec_Axe1; 
Glob_V1  = Glob_Vec_VAxe1Seg1;
Glob_V2  = Glob_Vec_VAxe1Seg2;
Seg1_U   = Angle.Joint.Vector1.LocCoord;
Seg2_U   = Angle.Joint.Vector2.LocCoord;
Seg1_V1  = Angle.V1Seg1.LocCoord;
Seg2_V1  = Angle.V1Seg2.LocCoord;


% auxiliar variables
dotV1U = dot(Seg1_V1, Seg1_U);
dotV2U = dot(Seg2_V1, Seg2_U);
Seg1_V1p = Seg1_V1 - dotV1U * Seg1_U;
Seg2_V1p = Seg2_V1 - dotV2U * Seg2_U;
ModV1p = norm(Seg1_V1p);
ModV2p = norm(Seg2_V1p); 

% calc constraints for rotation
Phi = mkCtrAngle(Glob_U, Glob_V1, Glob_V2, dotV1U, dotV2U, ModV1p, ModV2p, Alpha);
