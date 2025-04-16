function [qdot_t,qdot2_t] = calcVelAcc(q_t,Deltat)
% CALCVELACCEL calculates the velocity and acceleration of each
% generalized coordinate in q_t using cubic splines.
%
%   calcVelAccel(q_t,Deltat)

% sizes
[nSamples, nCoords] = size(q_t);

if nSamples <= 3
    warning('Not enough frames to estimate Vel and Acc');
    qdot_t = [];
    qdot2_t = [];
    return
end

tend = Deltat * (nSamples-1);
time = 0.0:Deltat:tend;

% piecewise cubic splines:
%   position:     A*x^3 +   B*x^2 +   C*x^1 +   D = 0,  [break(i), break(i+1)]
%   velocity:             3*A*x^2 + 2*B*x^1 +   C = 0,  [break(i), break(i+1)]
%   acceleration:                   6*A*x^1 + 2*B = 0,  [break(i), break(i+1)]

% piecewise polynomial for position
warning('off','all');
pp_pos = spline(time,q_t');
warning('on','all');

% y = q_t(:,1);
% xx = 0.0:Deltat/10:tend;
% yy = spline(time,y,xx);
% plot(x,y,'o',xx,pp_pos)
% plot(time,y,'o',xx,yy)
% coeficients of piecewise polynomial
A = pp_pos.coefs(:,1);
B = pp_pos.coefs(:,2);
C = pp_pos.coefs(:,3);
D = pp_pos.coefs(:,4);
zerosCol  = zeros(length(A), 1);

% piecewise polynomial for velocity
pp_vel = pp_pos;
pp_vel.coefs = [zerosCol, 3*A, 2*B, C];

% piecewise polynomial for acceleration
pp_accel = pp_pos;
pp_accel.coefs = [zerosCol, zerosCol, 6*A, 2*B];

% calc velocities and accelerations
qdot_t  = ppval(pp_vel,time)';
qdot2_t = ppval(pp_accel,time)';