function Phi = mkCtrFixMod(V,locV)

% MKCTR_FIXMOD creates a fixed modulus constraint
%
%   Phi = mkCtrFixMod(V,locV)
%   Inputs:
%     + V is the symbolic vector - symbolic array(3x1).
%     + locV are the local coordinates of the vector - double array(3x1)
%   Outputs:
%     + Phi is the symbolic fixed modulus vector constraint


SquaredModulus  = locV' * locV;
Phi = (V' * V) - SquaredModulus;
    