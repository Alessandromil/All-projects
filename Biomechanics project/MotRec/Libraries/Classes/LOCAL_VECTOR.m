classdef LOCAL_VECTOR < handle
    %LOCALVECTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Vector      % Handle to vector
        CoordName   % Name of local points SegName_PointNamex/y/z   char{3,1}
        LocCoord    % local coord of Vector [3x1]
    end
    
    methods
        function LV = LOCAL_VECTOR(Vector,SegName)
            LV.Vector = Vector;
            LV.CoordName = {[SegName,'_',Vector.Name,'x'];[SegName,'_',Vector.Name,'y'];[SegName,'_',Vector.Name,'z']};
        end
    end
    
end

