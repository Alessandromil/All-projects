function [varargout] = rot32(gamma,beta)
% This function admits one call
%
% A) It receives three Euler angles 32(zy) and
%    returns the rotation matrix 32(zy),
%
%       [R] = rot32(gamma,beta)
% NOTE:
%    Notice that rotation 321(zyx) in follower axis gives the
%    same rotation matrix than rotation 123(xyz) in fixed axis.
%    The rotation angle in X,Y and Z are the same for both cases
%    only the order is different. 


varargout{1} = rot_x(gamma) * rot_y(beta) ;
if nargout > 1
    error('Too many output expected')
end





