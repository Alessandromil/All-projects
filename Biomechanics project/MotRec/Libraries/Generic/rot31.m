function [varargout] = rot31(gamma,alpha)
% This function admits one call
%
% A) It receives three Euler angles 31(zx) and
%    returns the rotation matrix 31(zx),
%
%       [R] = rot31(gamma,alpha)
% NOTE:
%    Notice that rotation 321(zyx) in follower axis gives the
%    same rotation matrix than rotation 123(xyz) in fixed axis.
%    The rotation angle in X,Y and Z are the same for both cases
%    only the order is different. 


varargout{1} = rot_z(gamma) * rot_x(alpha);
if nargout > 1
    error('Too many output expected')
end





