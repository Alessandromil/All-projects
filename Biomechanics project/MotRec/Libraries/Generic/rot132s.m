function [varargout] = rot132s(varargin)
% This function admits two different calls
%
% A) It receives three Euler angles 132(xzy) and
%    returns the rotation matrix 132(xzy),
%
%       [R] = rot132s(alpha,gamma,beta)
%
% B) It receives the rotation matrix 132(xzy) and
%    returns three Euler angles 132(xzy), 
%
%       [alpha,gamma,beta] = rot132s(R)
%

if nargin == 3
    alpha = varargin{1};
    gamma = varargin{2};   
    beta  = varargin{3};
    varargout{1} = rot_x(alpha) * rot_z(gamma) * rot_y(beta);
    if nargout > 1
        error('Too many output expected')
    end
end

if nargin == 1
    R = varargin{1};
    epsilon = 1e-4;
    if abs( R(1,2) ) == 1 - epsilon  
        warning('gamma is close to 90 degrees. Close to singular position');
    end
	gamma = asin(-R(1,2));
	beta  = atan2(R(1,3), R(1,1));
	alpha = atan2(R(3,2), R(2,2));
    varargout{1} = alpha;
    varargout{2} = gamma;
    varargout{3} = beta;
    if nargout > 3
        error('Too many outputs expected')
    end
end
