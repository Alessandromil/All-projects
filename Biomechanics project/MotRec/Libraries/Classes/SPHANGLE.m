classdef SPHANGLE 
    %SPHANGLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name    % Name of Angle                     char
        Joint   % Name of Joint asociated to angle  char
        V1Seg1  % Handle to vector in Seg1          LOCAL_VECTOR
        V1Seg2  % Handle to vector in Seg2          LOCAL_VECTOR
        V2Seg1  % Handle to vector in Seg1          LOCAL_VECTOR
        V2Seg2  % Handle to vector in Seg2          LOCAL_VECTOR
        PosInq  % Pos of angle in q.                int
        RotSeq  % Rotation Sequence                 char
        
    end
    
    methods
        function SPHA = SPHANGLE(Name,Joint,V1Seg1,V1Seg2,V2Seg1,V2Seg2,RotSeq)
            SPHA.Name = Name;
            SPHA.Joint = Joint;
            SPHA.V1Seg1 = V1Seg1;
            SPHA.V1Seg2 = V1Seg2;
            SPHA.V2Seg1 = V2Seg1;
            SPHA.V2Seg2 = V2Seg2;
            SPHA.RotSeq = RotSeq;
        end
    end
    
end

