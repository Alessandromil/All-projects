function Var = mkSymbAngle(Name)

% MKSYMBANGLE creates a symbolic angle variable with name 'Name'
%
%   Var = mkSymbAngle(Name)
%   Inputs:
%     + Name is the name of the angle - string
%   Outputs:
%     + Var is a symbolic scalar variable(1x1) with the angle

Var = sym(Name,'real');
