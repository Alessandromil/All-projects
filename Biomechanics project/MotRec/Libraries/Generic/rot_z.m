function [R] = rot_z(gamma)
% Returns the rotation matriz in Z. Gamma can be
% numeric or symbolic. e.g.
%      [R] = rot_z(gamma) 


    R = [cos(gamma)  -sin(gamma)   0
         sin(gamma)   cos(gamma)   0
         0            0            1];
return