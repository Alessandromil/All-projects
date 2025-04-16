function [varargout] = rot232s(varargin)
% This function admits two different calls
%
% A) It receives three Euler angles 232(yzy) and
%    returns the rotation matrix 232(yzy),
%
%       [R] = rot232s(beta1,gamma,beta2)
%
% B) It receives the rotation matrix 232(yzy) and
%    returns three Euler angles 232(yzy), 
%
%       [beta1,gamma,beta2] = rot232s(R)
%

if nargin == 3
    beta1 = varargin{1};
    gamma = varargin{2};   
    beta2 = varargin{3};
    varargout{1} = rot_y(beta1) * rot_z(gamma) * rot_y(beta2);
    if nargout > 1
        error('Too many output expected')
    end
end

if nargin == 1
    error('This function version has not been implemented yet.')
end
