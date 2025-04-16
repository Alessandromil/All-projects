function [varargout] = rot21(beta,alpha)
% This function admits one call
%
% A) It receives three Euler angles 21(yx) and
%    returns the rotation matrix 21(zyx),
%
%       [R] = rot21(beta,alpha)
% NOTE:
%    Notice that rotation 321(zyx) in follower axis gives the
%    same rotation matrix than rotation 123(xyz) in fixed axis.
%    The rotation angle in X,Y and Z are the same for both cases
%    only the order is different. 


varargout{1} = rot_y(beta) * rot_x(alpha);
if nargout > 1
    error('Too many output expected')
end





