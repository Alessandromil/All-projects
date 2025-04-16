function [Deltaq, rankC] = spSolverQR2s(H,A,d,Tol)

% SPSOLVERQR2S solves a sparse linear system C * y = d using
% a sparse solver based on the QR decomposition in TWO STEPS (2s)
% See function QR for more details on QR decomposition.
%
%   Deltaq = spSolverQR2s(H,A,d,Tol)
%
%   Inputs:
%     + H is the Hessian matrix of the objective function
%     + A is Phiq'*Phiq
%     + d is the vector of independent terms
%     + Tol is the tolerance
%   Outputs:
%     + Deltaq is the solution vector
%     + rankC is the rank of the system matrix C


% The linear system of equations is:
%   [H  A][ Deltaq] = [-g] ->  C*y=d
%   [A  0][-Lambda]   [ b]
% Consider the following matrices
%   E=[H]  F=[A] d = [-g]
%     [A]    [0]     [ b]
% and that Lambda = -Lambda. The original variable can be obtained
% by changing the sign of the new Lambda
%
% Then 
%    E*Deltaq + F*Lambda = d [1]
%


% 1) Form matrices E, F, d -> E*deltaq + F*Lambda = d [1]   ----------------------------------
% sizes
nCoords = size(A,1);
% Definitions
E = [H; A];
F = [A; sparse(zeros(nCoords,nCoords))];
% d is an input of the function


% 2) QR decomposition of E:  E = QR ----------------------------------------------------------
% Note that orthogonal transformations are applied to [F d] producing
%  Q'*[F d] = [Q'*F  Q'*d] without computing Q.
%
% From [1] Q*R*Deltaq + F*Lambda = d
%            R*Deltaq + Q'*F*Lambda = Q'd
%            R*Deltaq + F0*Lambda = d0
[Fd,R]= qr(E,[F, d]); 
d0 = Fd(:,end);
F0 = Fd(:,1:end-1);


% 3) Analizing the structure of [1] we can get two linear system of equations  ---------------
%            [ R1 ] * deltaq + [ F1 ] * lambda = [ d1 ]
%            [  0 ]            [ F2 ]            [ d2 ]
%
%   R1*Deltaq + F1*lambda = d1 [2]
%               F2*lambda = d2 [3]

% Definition of matrices in linear system [2]
R1 = R(1:nCoords,:);
F1 = F0(1:nCoords,:);
d1 = d0(1:nCoords);

% Definition of matrices in linear system [3]
F2 = F0(nCoords+1:end,:);
d2 = d0(nCoords+1:end);

% Check if the problem is well formulated, i.e. all DoFs are effectively driven
diagR1    = full(diag(R1));
NoPivot   = find(abs(diagR1)<Tol);
nNoPivots = length(NoPivot);
if  nNoPivots > 0
    error(['The model is not properly guided. There are ',num2str(nNoPivots),' DoFs not guided']); 
end


% 4) QR decomposition of matrix F2 for linear system [3]  ------------------------------------
%    F2 = Q2*R2
%
%    Q2*R2*Lambda = d2        [3]
%    R2*Lambda = Q2'*d2 = d3  [3]
[d3, R2] = qr(F2, d2);


% 5) Analizing the structure of [3] ----------------------------------------------------------
%                                   [ R4| D4 ] * Lambda = [d4]
%     R2*Lambda = Q2'*d2 = d3 ->    [  0 |  0 ]            [ 0] 
%                                   
% If R2 is ordered as an upper triangular matrix then we can distinguish R3 and D3
% Note that after reordering the vector of Lambdas will change.
%
%    R2*Lambda = Q2'*d2 = d3 -> R2 * E2 * (E2' * Lambda) = d3
%
% with new variables:
%    Lambda2 = E2' * Lambda
%    R3      = R2  * E2
%
% Then
%    R3*Lambda2 = d3

% 5.a) Calculating permutation (column reordering) matrix E2 ***************
% size
[nRows_R2, nCols_R2] = size(R2);

% Initialize variables
pivotColIndex    = []; 
noPivotColIndex  = []; 
nPivots          =  0;

% Find Pivot column index and noPivot column index
for j = 1 : nCols_R2
    i = nPivots + 1;
    if(abs(R2(i,j)) < Tol)
        noPivotColIndex = [noPivotColIndex; j];
    else
        nPivots = nPivots + 1;
        pivotColIndex = [pivotColIndex; j];
    end
end

% permutation matrix. First: Pivots Cols. After noPivot Cols
E = speye(nCols_R2);
E2 = E(:,[pivotColIndex' noPivotColIndex']);


% 5.b) Calculate R3 = R2 * E2 ***********************************************
% After reordering we have R2*E2 = [ R3 | D3 ]
%                                  [  0 |  0 ]
R3 = R2 * E2; 


% 5.c) Extract vars R3, D3, b4 from R3 and b3 *******************************
R4 = R3(1:nPivots,1:nPivots);
%D4 = R3(nPivots+1:end,nPivots+1:end); % never used
d4 = d3(1:nPivots); 

% 5.d) Calculate lambda *****************************************************
Lambda2_free = zeros(nRows_R2-nPivots,1);
Lambda2_ind = R4\d4;
Lambda2 = [Lambda2_ind; Lambda2_free];
% undo the column(variable) reordering
Lambda = E2 * Lambda2;


% 6) From linear equation [2] Deltaq = R1\(d1-F1*Lambda) --------------------------------------------
Deltaq = full(R1\(d1-F1*Lambda));


% 7) Rank of the coefficient matrix  ----------------------------------------------------------------
% We suppose that E has full column rank otherwise and error is displayed.
% We know that E and F are L.I. then:
%     rank(C) = rank([E F]) = rank(E) + rank(F) = nCoords + rank(F) = nCoords + rank(R3)
rankC = nCoords + nPivots;

% 7) Check if the problem is under-guided -----------------------------------------------------------------
% if ReordMeth == 2 %#### 0
    nVars = nRows/2;
    nVarsUnder = 0;
    ValueVariable = [];
    for i=1:length(noPivotColIndex)
        if noPivotColIndex(i) <= nVars
            nVarsUnder = nVarsUnder + 1;
            ValueVariable = [ValueVariable;noPivotColIndex(i)];
        end
    end
    if nVarsUnder > 0
        str1 = sprintf(['    WARNING: There could be ',num2str(nVarsUnder),' variables under-guided !!!\n']);
        str2 = '    The underguide variables are: ';
        for i=1:nVarsUnder
            str2 = [str2,S.Reconstruction.Subject.q{ValueVariable(i)},' '];
        end
        str3 = [str1,str2];
        if(S.Display == 1 || S.Display == 2)
            disp(str3);
        end
        fprintf(S.ExperLogFileId, '%s\n', str3);
    end
% end