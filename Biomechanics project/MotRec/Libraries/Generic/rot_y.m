function [R] = rot_y(beta)
% Returns the rotation matriz in Y. Beta can be
% numeric or symbolic. e.g.
%      [R] = rot_y(beta) 

    R = [ cos(beta)   0   sin(beta)
          0           1   0
         -sin(beta)   0   cos(beta)];
return