function TrajSmark_filt = filtTraj(TrajSmark, SmoothPar)

% FILTTRAJ smooth the trajectory of one marker using method defined
% in cell SmoothPar
% Inputs:
%   TrajSmark: array(nFrames x 3) with x,y,z marker coordinates
%   SmoothPar cell{3x1}:
%      SmoothPar{1} string: name of the smoothing method
%      SmoothPar{2} array (nPar x 1): contain parameter values for the method
%      SmoothPar{3} sampling freq. of motion capture system
%
% Output:
%    TrajSmark_filt array(nFrames x 3) with smoothed marker trajectory


SmoothMethod = SmoothPar{1}; % smoothing method. string
MethodPar    = SmoothPar{2}; % smoothing method parameters. array(nPar x 1)
sampleFreq   = SmoothPar{3}; % sampling freq. of the motion capture system



if strcmpi(SmoothMethod,'spline')
    % Spline smoothing
    p = MethodPar; 
    TrajSmark_filt = filtSpline(TrajSmark, p);

elseif strcmpi(SmoothMethod,'butter')
    % Butterworth filtering 
    cutFreq = MethodPar; 
    TrajSmark_filt = filtButter(TrajSmark, sampleFreq, cutFreq);
    
%elseif strcmpi(SmoothMethod,'ggpsa')  -> Sergio
%elseif strcmpi(SmoothMethod,'fft')    -> Sergio
%elseif strcmpi(SmoothMethod,'winter') -> Francois Fraisse

end

