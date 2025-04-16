function [R] = rot_x(alpha)
% Returns the rotation matriz in X. Alpha can be
% numeric or symbolic. e.g.
%      [R] = rot_x(alpha) 

    R = [1   0            0
         0   cos(alpha)  -sin(alpha)
         0   sin(alpha)   cos(alpha)];
return