function Sym = mkSymbolicXYZ(Name)
% MKSYMBOLIC returns a symbolic vector (3x1) with the 3D coordinates
% of 'Name'
%
%   Sym = mkSymbolicXYZ(Name)
%   Inputs:
%     + Name is the name of the variable - string
%   Outputs:
%     + Sym is a symbolic array (3x1) with the x,y and z coordinates
%       of the symbolic variable
%
%   Example:
%      Vector1 = mkSymbolicXYZ('V1');
%      where
%       Vector1(1) = V1x (symbolic variable)
%       Vector1(2) = V1y (symbolic variable)
%       Vector1(3) = V1z (symbolic variable)

x = sym([Name,'x'], 'real');
y = sym([Name,'y'], 'real');
z = sym([Name,'z'], 'real');

Sym = [x; y; z];