function [varargout] = rot312s(varargin)
% This function admits two different calls
%
% A) It receives three Euler angles 312(zxy) and
%    returns the rotation matrix 312(zxy),
%
%       [R] = rot231s(gamma,alpha,beta)
%
% B) It receives the rotation matrix 312(zxy) and
%    returns three Euler angles 312(zxy), 
%
%       [gamma,alpha,beta] = rot312s(R)
%

if nargin == 3
    gamma = varargin{1};   
    alpha = varargin{2};
    beta  = varargin{3};
    varargout{1} = rot_z(gamma) * rot_x(alpha) * rot_y(beta);
    if nargout > 1
        error('Too many output expected')
    end
end

if nargin == 1
    R = varargin{1};
    epsilon = 1e-4;
    if abs( R(3,2) ) == 1 - epsilon  
        warning('alpha is close to 90 degrees. Close to singular position');
    end
	alpha = asin(R(3,2));
	beta  = atan2(-R(3,1), R(3,3));
	gamma = atan2(-R(1,2), R(2,2));
    varargout{1} = gamma;
    varargout{2} = alpha;
    varargout{3} = beta;
    if nargout > 3
        error('Too many outputs expected')
    end
end
