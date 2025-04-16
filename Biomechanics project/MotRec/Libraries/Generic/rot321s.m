function [varargout] = rot321s(varargin)
% This function admits two different calls
%
% A) It receives three Euler angles 321(zyx) and
%    returns the rotation matrix 321(zyx),
%
%       [R] = rot321s(gamma,beta,alpha)
%
% B) It receives the rotation matrix 321(zyx) and
%    returns three Euler angles 321(zyx), 
%
%       [gamma,beta,alpha] = rot321s(R)
%
% NOTE:
%    Notice that rotation 321(zyx) in follower axis gives the
%    same rotation matrix than rotation 123(xyz) in fixed axis.
%    The rotation angle in X,Y and Z are the same for both cases
%    only the order is different. 

if nargin == 3
    gamma = varargin{1};
    beta  = varargin{2};
    alpha = varargin{3};   
    varargout{1} = rot_z(gamma) * rot_y(beta) * rot_x(alpha);
    if nargout > 1
        error('Too many output expected')
    end
end

if nargin == 1
    R = varargin{1};
    epsilon = 1e-4;
    if abs( R(3,1) ) == 1 - epsilon  
        warning('beta is close to 90 degrees. Close to singular position');
    end
	beta  = asin(-R(3,1));
	alpha = atan2(R(3,2), R(3,3));
	gamma = atan2(R(2,1), R(1,1));
    varargout{1} = gamma;
    varargout{2} = beta;
    varargout{3} = alpha;
    if nargout > 3
        error('Too many outputs expected')
    end
end


