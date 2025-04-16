function Phi = mkCtrSPHAngle(Angle)

% MKCTRSPHANGLE creates the relative angle constraints for a spheric joint
%
%   Phi = mkCtrSPHAngle(SPHAngle)
%   Inputs:
%     + SPHAngle is a class with other to do the constraint
%   Outputs:
%     + Phi is the symbolic Universal joint constraint

% set variables
Sequence = Angle.RotSeq;
j_prev = Angle.V1Seg1.Vector.Name;
k_prev = Angle.V2Seg1.Vector.Name;
j_next = Angle.V1Seg2.Vector.Name;
k_next = Angle.V2Seg2.Vector.Name;

% ### Por el momento lo hago así luego cuando tengamos el modelo de ESI
% vemos si tenemos giros que no se correspondan con los ejes principales
alpha = Angle.Name1;  % relative Euler angle between Body1 and Body2 in x-axis
beta  = Angle.Name2;  % relative Euler angle between Body1 and Body2 in y-axis
gamma = Angle.Name3;  % relative Euler angle between Body1 and Body2 in z-axis
% if strcmp(Sequence,'123')
%     alpha = Angle.Name1;  % relative Euler angle between Body1 and Body2 in x-axis
%     beta  = Angle.Name2;  % relative Euler angle between Body1 and Body2 in y-axis
%     gamma = Angle.Name3;  % relative Euler angle between Body1 and Body2 in z-axis
% elseif strcmp(Sequence,'321')
%     gamma = Angle.Name1;  % relative Euler angle between Body1 and Body2 in z-axis
%     beta  = Angle.Name2;  % relative Euler angle between Body1 and Body2 in y-axis
%     alpha = Angle.Name3;  % relative Euler angle between Body1 and Body2 in x-axis
% elseif strcmp(Sequence,'231')
%     beta  = Angle.Name1;  % relative Euler angle between Body1 and Body2 in y-axis
%     gamma = Angle.Name2;  % relative Euler angle between Body1 and Body2 in z-axis
%     alpha = Angle.Name3;  % relative Euler angle between Body1 and Body2 in x-axis
% elseif strcmp(Sequence,'232')
%     beta1 = Angle.Name1;  % relative Euler angle between Body1 and Body2 in y-axis
%     gamma = Angle.Name2;  % relative Euler angle between Body1 and Body2 in z-axis
%     beta2 = Angle.Name3;  % relative Euler angle between Body1 and Body2 in y-axis
% elseif strcmp(Sequence,'132')
%     alpha = Angle.Name1;  % relative Euler angle between Body1 and Body2 in x-axis
%     gamma = Angle.Name2;  % relative Euler angle between Body1 and Body2 in z-axis
%     beta  = Angle.Name3;  % relative Euler angle between Body1 and Body2 in y-axis
% elseif strcmp(Sequence,'312')
%     gamma = Angle.Name1;  % relative Euler angle between Body1 and Body2 in z-axis
%     alpha = Angle.Name2;  % relative Euler angle between Body1 and Body2 in x-axis
%     beta  = Angle.Name3;  % relative Euler angle between Body1 and Body2 in y-axis
% else
%     error(sprintf(['\n------------------------------------------------\n'...
%                    ' Rotation sequence not defined or wrong format !!'...
%                    '\n------------------------------------------------\n']));
% end

% Check each input. If it is a char variable create the corresponding symbolic variable
if ischar(j_prev), j_prev = mkSymbolicXYZ(j_prev);  end
if ischar(k_prev), k_prev = mkSymbolicXYZ(k_prev);  end
if ischar(j_next), j_next = mkSymbolicXYZ(j_next);  end
if ischar(k_next), k_next = mkSymbolicXYZ(k_next);  end
if ischar(alpha),  alpha  = mkSymbAngle(alpha); end
if ischar(beta),   beta   = mkSymbAngle(beta);  end
if ischar(gamma),  gamma  = mkSymbAngle(gamma); end
 

% relative rotation matrix between body 1 (previous) and body 2 (next)
if strcmp(Sequence,'123')
    Rot = rot123s(alpha, beta, gamma);
elseif strcmp(Sequence,'321')
    Rot = rot321s(gamma, beta, alpha);
elseif strcmp(Sequence,'231')
    Rot = rot231s(beta, gamma, alpha);
elseif strcmp(Sequence,'232')
    Rot = rot232s(beta1, gamma, beta2);
elseif strcmp(Sequence,'132')
    Rot = rot132s(alpha, gamma, beta);    
elseif strcmp(Sequence,'312')
    Rot = rot312s(gamma, alpha, beta);    
else
    error(sprintf(['\n------------------------------------------------\n'...
                   ' Rotation sequence not defined or wrong format !!'...
                   '\n------------------------------------------------\n']));
end

% rotation matrix between Ground (or World) and body 1 (previous)
i_prev = cross(j_prev, k_prev);
Glob_R_Prev = [i_prev, j_prev, k_prev];

% rotation matrix between Ground (or World) and body 2 (next)
i_next = cross(j_next, k_next);
Glob_R_Next = [i_next, j_next, k_next];

% rotation matrix between body 1 (previous) and body 2 (next)
Prev_R_Next = Glob_R_Prev' * Glob_R_Next;

% constraint equations 
k = 0;
Phi = sym('Phi','real');
% for j = 1:2
for j = 2:3
    for i = 1:3
        k = k + 1;
        Phi(k,1) = Rot(i,j) - Prev_R_Next(i,j);
    end
end



