function Phi = mkCtrUNIAngle(Angle)

% MKCTRANGLE creates the relative angle constraints for a Universal joint
%
%   Phi = mkCtrAngle(Angle)
%   Inputs:

%   Outputs:
%     + Phi are the symbolic Revolute joint constraints

% set variables
Axe1  = Angle.Joint.Vector1.Vector.Name; % First Axis of the Universal joint & reference vector in Seg2 It MUST be a UNITARY Vector
Axe2  = Angle.Joint.Vector2.Vector.Name; % Second Axis of the Universal joint & reference vector in Seg1 It MUST be a UNITARY Vector
VAxe1Seg1 = Angle.V1Seg1.Vector.Name; % Reference vector in Seg1. It MUST be a UNITARY Vector
VAxe2Seg2 = Angle.V1Seg2.Vector.Name; % Reference vector in Seg2. It MUST be a UNITARY Vector
% VAxe2Seg2 = Angle.V2Seg2.Vector.Name; % Reference vector in Seg2. It MUST be a UNITARY Vector
% VAxe2Seg1 = Angle.V2Seg1.Vector.Name; % Reference vector in Seg1. It MUST be a UNITARY Vector
Axe1_ang_Axe2 = Angle.Joint.CAngle; % Angle between Axe1 and Axe2. It must be constant (Ujoint constraint)
Alpha = Angle.Name1; % First relative angle of the UNIJoint
Beta  = Angle.Name2; % Second relative angle of the UNIJoint


% Check each input. If it is a char variable create the corresponding symbolic variable
if ischar(Axe1), Glob_Vec_Axe1 = mkSymbolicXYZ(Axe1);  end
if ischar(Axe2), Glob_Vec_Axe2 = mkSymbolicXYZ(Axe2);  end
if ischar(VAxe1Seg1), Glob_Vec_VAxe1Seg1 = mkSymbolicXYZ(VAxe1Seg1);  end
if ischar(VAxe2Seg2), Glob_Vec_VAxe2Seg2 = mkSymbolicXYZ(VAxe2Seg2);  end
% if ischar(VAxe2Seg2), Glob_Vec_VAxe2Seg2 = mkSymbolicXYZ(VAxe2Seg2);  end
% if ischar(VAxe2Seg1), Glob_Vec_VAxe2Seg1 = mkSymbolicXYZ(VAxe2Seg1);  end
if ischar(Alpha), Alpha = mkSymbAngle(Alpha);  end
if ischar(Beta),  Beta = mkSymbAngle(Beta);  end

% -----------------------------------------------
% FIRST Rotation
% -----------------------------------------------

Glob_U   = Glob_Vec_Axe1; 
Glob_V1  = Glob_Vec_VAxe1Seg1;
Glob_V2  = Glob_Vec_Axe2;
Seg1_U   = Angle.Joint.Vector1.LocCoord;
Seg1_V1  = Angle.V1Seg1.LocCoord;
U_ang_V2 = Axe1_ang_Axe2;

% auxiliar variables
dotV1U = dot(Seg1_V1, Seg1_U);
dotV2U = cos(U_ang_V2);
B1_V1p = Seg1_V1 - dotV1U * Seg1_U;
%B2_V2p = B2_V2 - dotV2U * B2_U; % It is no possible to get this value in local coordinates
ModV1p = norm(B1_V1p);
ModV2p = sin(U_ang_V2); 

% calc constraints for 1st rotation
phi_alpha = mkCtrAngle(Glob_U, Glob_V1, Glob_V2, dotV1U, dotV2U, ModV1p, ModV2p, Alpha);

% -----------------------------------------------
% SECOND Rotation
% -----------------------------------------------
Glob_U   = Glob_Vec_Axe2; 
Glob_V1  = Glob_Vec_VAxe2Seg2;
Glob_V2  = Glob_Vec_Axe1;
Seg2_U     = Angle.Joint.Vector2.LocCoord;
Seg2_V2    = Angle.V1Seg2.LocCoord;
V1_ang_U = Axe1_ang_Axe2;

% auxiliar variables
dotV1U = cos(V1_ang_U);
dotV2U = dot(Seg2_U, Seg2_V2);
%B1_V1p = B1_V1 - dotV1U * B1_U; % It is no possible to get this value in local coordinates
B2_V2p = Seg2_V2 - dotV2U * Seg2_U; 
ModV1p = sin(V1_ang_U); 
ModV2p = norm(B2_V2p); 

% calc constraints for 1st rotation
phi_beta = mkCtrAngle(Glob_U, Glob_V1, Glob_V2, dotV1U, dotV2U, ModV1p, ModV2p, Beta);

% output
Phi = [phi_alpha; phi_beta];

