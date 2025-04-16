function q0 = velapprox(q_t, i)

% VELAPPROX give the next initial approximation using the 
% velocity approximation method
%
%  q0 = velapprox(q_t, i)
%   Inputs:
%     + q_t is a double array(nFrames x nVars) with the value of each
%       generalized coordinate at each sample time. For the sample time
%       i, q(i,:) containts the values of the generalized coordinates.
%       The order of the variables is the same as in q.
%     + i is the number of the current Newton-Raphson step.  
%   Outputs:
%     + q0 is the initial approximation of the vector of generalized 
%       coordinates obtained from the velocity approximation


if i == 1
%     disp('First and second step without velocity approximation');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FALTA LOG
    q0 = q_t(i,:)';
else
    q0 = 2*q_t(i,:)' - q_t(i-1,:)';
end
