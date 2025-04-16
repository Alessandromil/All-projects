  function Phi = mkCtrLinComb(U, locU, SegRef,InvR)

% MKCTR_LINCOMB creates a linear combination constraint. This function
% express a point or a vector as a linear combination of a reference
% and creates the corresponding constraint
% 
%   Phi = mkCtrLinComb(U, locU, SymAxes)
%   Inputs:
%     + U is a cell (nElements x 1) with the symbolic elements (points or vectors)
%       that belong to a body
%     + locU is a cell (nElements x 1) with the local coordinates of the elements
%       (points or vectors) that belong to a body 
%     + SegRef is a struct (1x1) that contain:
%         Diri: strruc (1x1) that contain:
%               Global : sym value (3x1)
%               Local  : sym or double value (3x1)
%     + InvR is the inverse of the projection of vector U on each axis of the reference.
%       R = [SegRef.Dir1.Local, SegRef.Dir2.Local, SegRef.Dir3.Local]
%   Outputs:
%     + Phi are the symbolic linear combination constraints


% coeficients
coef = InvR * locU;

% contraint equations
Phi = U - coef(1) * SegRef.Dir1.Global - coef(2) * SegRef.Dir2.Global - coef(3) * SegRef.Dir3.Global; % 3 constraint equations