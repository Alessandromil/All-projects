function Phi = mkCtrUjoint(varargin)

% MKCTRUJOINT creates a Universal joint constraint. The inputs arguments
% are different depending on the joint configuration.
% 
%   Phi = mkctr_Ujoint(varargin)
%
%   If the two vectors that define a Universal joint are perpendicular:
%   Inputs:
%     + Vec_i is the LOCAL_VECTOR that defines the 1st axis of the Ujoint -
%       in Seg1
%     + Vec_j is the LOCAL_VECTOR that defines the 2nd axis of the Ujoint -
%       in Seg2
%
%   If the two vectors that define a Universal joint are NOT perpendicular:
%   Inputs (OPTION 1):
%     + Vec_i is the LOCAL_VECTOR that defines the 1st axis of the Ujoint -
%       in Seg1
%     + Vec_j is the LOCAL_VECTOR that defines the 2nd axis of the Ujoint -
%       in Seg2
%     + Ang is the angle between both axis of the Ujoint - scalar
%   Inputs (OPTION 2): No adapted to class
%     + Vec_i is the symbolic vector that defines the 1st axis of the Ujoint -
%       symbolic array(3x1)
%     + Vec_j is the symbolic vector that defines the 2nd axis of the Ujoint -
%       symbolic array(3x1)
%     + Loc_i are the local coordinates of the 1st axis of the Ujoint -
%       double array (3x1). The local coordinates of the 1st and 2nd axes must
%       be given in the same coordinate system
%     + Loc_j are the local coordinates of the 2nd axis of the Ujoint -
%       double array (3x1). The local coordinates of the 1st and 2nd axes must
%       be given in the same coordinate system
%
%   Outputs:
%     + Phi is the symbolic Universal joint constraint



if nargin == 2
    % both vectors are perpendicular
    Vec_i = varargin{1}.Vector.Name;
    Vec_j = varargin{2}.Vector.Name;
    % check if Vec_i and Vec_j are a symbolic array or a char
    if ischar(Vec_i), Vec_i = mkSymbolicXYZ(Vec_i); end
    if ischar(Vec_j), Vec_j = mkSymbolicXYZ(Vec_j); end
    % build constraint
    Phi = dot(Vec_i, Vec_j);
  
elseif nargin == 3
    % The vectors are NOT perpendicular. The angle is given
    Vec_i = varargin{1}.Name;
    Vec_j = varargin{2}.Name;
    Ang   = varargin{3}; % [rad]
    % check if Vec_i and Vec_j are a symbolic array or a char
    if ischar(Vec_i), Vec_i = mkSymbolicXYZ(Vec_i); end
    if ischar(Vec_j), Vec_j = mkSymbolicXYZ(Vec_j); end
    % build constraint
    Phi = dot(Vec_i, Vec_j) - cos(Ang);
    % warnings
    tmp = char(Vec_i(1)); AxisName_i = tmp(1:end-1);
    tmp = char(Vec_j(1)); AxisName_j = tmp(1:end-1);
    str = sprintf(['\n---------------------------------------------------------------------------------------------------------\n', ...
                   ' Universal joint defined with axes "',AxisName_i,'" and "',AxisName_j,'": Both vectors MUST have modulus 1.0\n', ...
                   '---------------------------------------------------------------------------------------------------------\n']);
    fprintf(MODEL_LOGFILE_ID, '%s\n', str);
    if INFO_TYPE == 1
        disp(str);
    end
end

% elseif nargin == 4
%     % The vectors are NOT perpendicular.
%     % local coord. of both vectors are given in a common frame.
%     Vec_i = varargin{1};
%     Vec_j = varargin{2};
%     Loc_i = varargin{3};
%     Loc_j = varargin{4};
%     % check if Vec_i and Vec_j are a symbolic array or a char
%     if ischar(Vec_i), Vec_i = mkvector(Vec_i); end
%     if ischar(Vec_j), Vec_j = mkvector(Vec_j); end
%     % build constraint
%     Phi = dot(Vec_i, Vec_j) - norm(Loc_i)*norm(Loc_j)*dot(Loc_i, Loc_j);    
% end
