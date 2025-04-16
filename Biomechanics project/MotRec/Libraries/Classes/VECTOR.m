classdef VECTOR < handle
    %VECTOR Summary of this class goes here
        
    properties
        Name        % Name of the vector.           char
        CoordName   % Name of the vector coord      char{3,1}
        GlobalCoord % Global coor of vector.        double[3x1]
        PosInq      % Pos of vector in q.           int
        Fixed = 0;  % Is in fixed segment Fixed==1
    end
    
    methods
        function V = VECTOR(Name)
            V.Name = Name;
            V.CoordName = {[Name,'x'];[Name,'y'];[Name,'z']};
        end
        
    end
    
end

