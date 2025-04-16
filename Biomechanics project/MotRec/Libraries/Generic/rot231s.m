function [varargout] = rot231s(varargin)
% This function admits two different calls
%
% A) It receives three Euler angles 231(yzx) and
%    returns the rotation matrix 231(yzx),
%
%       [R] = rot231s(beta,gamma,alpha)
%
% B) It receives the rotation matrix 231(yzx) and
%    returns three Euler angles 231(yzx), 
%
%       [beta,gamma,alpha] = rot231s(R)
%

if nargin == 3
    beta  = varargin{1};
    gamma = varargin{2};   
    alpha = varargin{3};
    varargout{1} = rot_y(beta) * rot_z(gamma) * rot_x(alpha);
    if nargout > 1
        error('Too many output expected')
    end
end

if nargin == 1
    R = varargin{1};
    epsilon = 1e-4;
    if abs( R(2,1) ) == 1 - epsilon  
        warning('gamma is close to 90 degrees. Close to singular position');
    end
	gamma = asin(R(2,1));
	beta  = atan2(-R(3,1), R(1,1));
	alpha = atan2(-R(2,3), R(2,2));
    varargout{1} = beta;
    varargout{2} = gamma;
    varargout{3} = alpha;
    if nargout > 3
        error('Too many outputs expected')
    end
end
