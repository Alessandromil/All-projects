function Phi = mkCtrFixAng(varargin)

% MKCTR_FIXANG creates a fixed angle constraint between two vectors
%
%   Phi = mkCtrFixAng(varargin)
%   Inputs:
%     + V1 is the symbolic vector 1 - symbolic array(3x1).
%     + V2 is the symbolic vector 2 - symbolic array(3x1).
%     + locV1 are the local coordinates of the vector 1 - double array(3x1) or symbolic array (3x1)
%     + locV2 are the local coordinates of the vector 2 - double array(3x1) or symbolic array (3x1)
%     + epsilon (OPTIONAL ARGUMENT) tolerance to consider cos(angle) = 0 
%   Outputs:
%     + Phi is the symbolic fixed angle constraint


% check number of input/output arguments
if nargin == 4
    epsilon = 1.0e-12;
elseif nargin == 5
    epsilon = varargin{5};
else    
    error('Number of arguments is incorrect')
end
if nargout > 1
    error('Too many output expected')
end

% set variables
V1 = varargin{1};
V2 = varargin{2};
locV1 = varargin{3};
locV2 = varargin{4};

% Check if Loc are symbolic
if isobject(locV1) || isobject(locV2)
    % vector lengths
    LengthV1 = normSym(locV1);
    LengthV2 = normSym(locV2); 
   
    % constant relative angle
    alpha  = acos( dot(locV1,locV2) / (LengthV1*LengthV2) );
    
    % constraint equation
    Phi = dot(V1,V2) - LengthV1 * LengthV2 * cos(alpha);
else    
    % vector lengths
    LengthV1 = norm(locV1);
    LengthV2 = norm(locV2);
    
    % constant relative angle
    alpha  = acos( dot(locV1,locV2) / (LengthV1*LengthV2) );
    
    % contraint equation
    if abs(cos(alpha)) < epsilon        % when alpha ~ 90 -> cos(alpha) ~ 0
        Phi = dot(V1,V2);               % alpha is assumed to be ZERO
    else
        Phi = dot(V1,V2) - LengthV1 * LengthV2 * cos(alpha);
    end
end

