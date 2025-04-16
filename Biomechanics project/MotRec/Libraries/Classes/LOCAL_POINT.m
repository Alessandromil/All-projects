classdef LOCAL_POINT < handle
    %LOCALPOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name        % the name of the point                         char
        Point       % Handle to                                     POINT
        CoordName   % Name of local points SegName_PointNamex/y/z   char{3,1}
        LocCoord    % local coord of Point.                         double[3x1]
        Origin = 0; % If local point is origin of Segment Origin=1. double[1x1]
    end
    
    methods
        function LP = LOCAL_POINT(Point,SegName)
            LP.Point = Point;
            LP.Name = Point.Name;
            LP.CoordName = {[SegName,'_',Point.Name,'x'];[SegName,'_',Point.Name,'y'];[SegName,'_',Point.Name,'z']};
        end
    end
    
end

