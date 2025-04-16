function Phi = mkCtrRevJoint(Vec1,Vec2)

Vec_i = Vec1.Vector.Name;
Vec_j = Vec2.Vector.Name;

% check if Vec_i and Vec_j are a symbolic array or a char
if ischar(Vec_i), Vec_i = mkSymbolicXYZ(Vec_i); end
if ischar(Vec_j), Vec_j = mkSymbolicXYZ(Vec_j); end

% build constraint
Phi = Vec_i - Vec_j;

end