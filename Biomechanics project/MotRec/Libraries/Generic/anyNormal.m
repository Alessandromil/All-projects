function VecN = anyNormal(Vec)

% ANYNORMAL returns any normal to the input vector.
%
%   VecN = anyNormal(Vec)
%   Inputs:
%     + Vec is 3D vector - double array(3x1) or (1x3)
%   Outputs:
%     + VecN is 3D vector perpendicular to the input
%       vector Vec - double array(3x1)


[val, index] = max(abs(Vec));

if index == 1
    VecN = [ Vec(2); -Vec(1);     0.0];
elseif index == 2
    VecN = [    0.0;  Vec(3); -Vec(2)];
elseif index == 3
    VecN = [-Vec(3);     0.0;  Vec(1)];
end

% normalize
VecN = VecN / norm(VecN);
