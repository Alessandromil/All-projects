classdef UNIJOINT < JOINT
    %UNIJOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Point1  % Handle to Local point of Seg1  LOCAL_POINT
        Point2  % Handle to Local point of Seg2  LOCAL_POINT
        Vector1   % Handle to vector in Seg1         LOCAL_VECTOR
        Vector2   % Handle to vector in Seg2         LOCAL_VECTOR
        CAngle    % Constant angle between V1 and V2 double
    end
    
    methods
        function UNIJ = UNIJOINT(Name,Type,Seg1,Seg2,Point1,Point2,Vector1,Vector2,CAngle)
            UNIJ = UNIJ@JOINT(Name,Type,Seg1,Seg2); % Call JOINT constructor
            UNIJ.Point1 = Point1;
            UNIJ.Point2 = Point2;
            UNIJ.Vector1 = Vector1;
            UNIJ.Vector2 = Vector2;
            UNIJ.CAngle = str2double(CAngle)*pi/180;
        end
    end
    
end

