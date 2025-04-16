function phi = mkCtrAngle(Glob_U, Glob_V1, Glob_V2, dotV1U, dotV2U, ModV1p, ModV2p, Theta)

% MKCTRANGLE creates the relative angle constraints for an angle
%
%   Phi = mkCtrAngle(Angle)
%   Inputs:

%   Outputs:
%     + Phi are the symbolic Revolute joint constraints


% auxiliar variables
crossV1V2 = cross(Glob_V1,Glob_V2);
crossV1U  = cross(Glob_V1,Glob_U);
crossUV2  = cross(Glob_U,Glob_V2);

% contraint equations
phi = sym('phi','real');
% dot constraint (x1)
phi(1,1) = dot(Glob_V1, Glob_V2) - dotV1U * dotV2U - ModV1p * ModV2p * cos(Theta);
% cross constraints (x3)
phi(2,1) = crossV1V2(1) - dotV2U*crossV1U(1) - dotV1U*crossUV2(1) - Glob_U(1) * ModV1p * ModV2p * sin(Theta); % ONLY ONE EQUATION IS NEEDED.
