classdef SPHJOINT < JOINT
    %SPHJOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Point1  % Handle to Local point of Seg1  LOCAL_POINT
        Point2  % Handle to Local point of Seg2  LOCAL_POINT
    end
    
    methods
        function SPHJ = SPHJOINT(Name,Type,Seg1,Seg2,Point1,Point2)
            SPHJ = SPHJ@JOINT(Name,Type,Seg1,Seg2); % Call JOINT constructor
            SPHJ.Point1 = Point1;
            SPHJ.Point2 = Point2;
        end
    end
    
end

