classdef ANGLE < handle
    %ANGLE In the Angle Class is the necesary information to make the
    %constrains depending the type of Joint.
    % For SPHJoint we need 3 Names and 4 vectors
    % For UNIJoint we need 2 Names and 4 vectors
    % For REVJoint we need 1 Name and 2 vectors
    
    properties
        Name1   % Name of Angle                     char
        Name2   % Name of Angle                     char
        Name3   % Name of Angle                     char
        a1      % Value of a1                       double
        a2      % Value of a2                       double
        a3      % Value of a3                       double
        Joint   % Handle to Joint                   JOINT
        V1Seg1  % Handle to vector in Seg1          LOCAL_VECTOR
        V1Seg2  % Handle to vector in Seg2          LOCAL_VECTOR
        V2Seg1  % Handle to vector in Seg1          LOCAL_VECTOR
        V2Seg2  % Handle to vector in Seg2          LOCAL_VECTOR
        PosInq  % Pos of angle in q.                int/int[]
        RotSeq  % Rotation Sequence                 char 
        
    end
    
    methods
        function A = ANGLE(Name1,Name2,Name3,Joint,V1Seg1,V1Seg2,V2Seg1,V2Seg2,RotSeq)
            A.Name1 = Name1;
            A.Name2 = Name2;
            A.Name3 = Name3;
            A.Joint = Joint;
            A.V1Seg1 = V1Seg1;
            A.V1Seg2 = V1Seg2;
            A.V2Seg1 = V2Seg1;
            A.V2Seg2 = V2Seg2;
            A.RotSeq = RotSeq;
        end
    end
    
end

