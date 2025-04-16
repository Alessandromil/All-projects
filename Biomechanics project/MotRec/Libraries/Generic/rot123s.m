function [varargout] = rot123s(varargin)

% This function admits two different calls
%
% A) It receives three Euler angles 123(xyz) and
%    returns the rotation matrix 123(xyz),
%
%       [R] = rot123s(alpha,beta,gamma)
%
% B) It receives the rotation matrix 123(xyz) and
%    returns three Euler angles 123(xyz), 
%
%       [alpha,beta,gamma] = rot123s(R)
%
% NOTE:
%    Notice that rotation 123(xyz) in follower axis gives the
%    same rotation matrix than rotation 321(zyx) in fixed axis.
%    The rotation angle in X,Y and Z are the same for both cases
%    only the order is different. 

if nargin == 3
    alpha = varargin{1};
    beta  = varargin{2};
    gamma = varargin{3};   
    varargout{1} = rot_x(alpha) * rot_y(beta) * rot_z(gamma);
    if nargout > 1
        error('Too many output expected')
    end
end

if nargin == 2
    R = varargin{1};
    PreviousGamma = varargin{2};
    epsilon = 1e-4;
    if abs( R(1,3) ) == 1 - epsilon  
        warning('beta is close to 90 degrees. Close to singular position');
    end
	beta  = asin(R(1,3));
	alpha = atan2(-R(2,3), R(3,3));
%     gamma = asin(-R(1,2)/cos(beta));
%     cosgamma = (R(1,1)/cos(beta));
%     if cosgamma < 0
%         if gamma > 0
%             gamma = pi - gamma;
%         elseif gamma < 0
%             gamma = -pi - gamma;
%         end
%     end
    CorrectAngle = 0;
	Gamma = atan2(-R(1,2), R(1,1));
    while CorrectAngle == 0
        GammaPlus = Gamma + 2*pi;
        GammaMinus = Gamma - 2*pi;
        DifGamma(1) = Gamma - PreviousGamma; DifGamma(2) = GammaPlus - PreviousGamma; DifGamma(3) = GammaMinus - PreviousGamma;
        if abs(DifGamma(2)) < abs(DifGamma(1)) && abs(DifGamma(2)) < abs(DifGamma(3))
            Gamma = GammaPlus;
        elseif abs(DifGamma(3)) < abs(DifGamma(1)) && abs(DifGamma(3)) < abs(DifGamma(2))
            Gamma = GammaMinus;
        else
            CorrectAngle = 1;
        end
    end
% gamma = atan(-R(1,2)/R(1,1));
%     if gamma <= -pi
%         gamma = -gamma;
%     end
    varargout{1} = alpha;
    varargout{2} = beta;
    varargout{3} = Gamma;
    if nargout > 3
        error('Too many outputs expected')
    end
end