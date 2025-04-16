classdef REVJOINT < JOINT
    %REVJOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Point1       % Handle to Local point of Seg1  LOCAL_POINT
        Point2       % Handle to Local point of Seg2  LOCAL_POINT
        Vector1      % Handle to vector in Seg1         LOCAL_VECTOR
        Vector2      % Handle to vector in Seg2         LOCAL_VECTOR
    end
    
    methods
        function REVJ = REVJOINT(Name,Type,Seg1,Seg2,Point1,Point2,Vector1,Vector2)
            REVJ = REVJ@JOINT(Name,Type,Seg1,Seg2); % Call JOINT constructor
            REVJ.Point1 = Point1;
            REVJ.Point2 = Point2;
            REVJ.Vector1 = Vector1;
            REVJ.Vector2 = Vector2;
        end
    end
    
end

