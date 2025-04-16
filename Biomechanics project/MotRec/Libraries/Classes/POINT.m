classdef POINT < handle
    %POINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        Name                % Name of the point.                                char
        CoordName           % Name of the vector coord                          char{3,1}
        GlobalCoord         % Global coor of point.                             double[3x1]
        PosInq              % Pos of point in q.                                int
        PosInz              % Pos of point in z.                                int
        Fixed = 0;          % Is in fixed segment                               Fixed = 1
        MeasuredCoord       % Global coord of markers                           double[3xNLandmarks]
        Postures            % Is the vector of the postures
                            % Posture.Name = The name of the posture            Char
                            % Posture.Glob = Global coord of posture            double[3x1]
        Wm                  % Weighting factor of the point                     double
        MarkerTrajectory    % Marker trajectory                                 double[NFramex3]
        OptPose = 0;        % It will be 1 when performing the optimization 
        
    end
    
    methods
        function P = POINT(Name)
            P.Name = Name;
            P.CoordName = {[Name,'x'];[Name,'y'];[Name,'z']};
        end

    end
    
end

