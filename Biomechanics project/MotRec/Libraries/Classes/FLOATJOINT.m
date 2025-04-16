classdef FLOATJOINT < JOINT
    %SPHJOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
      TrsSensor  % Handle to the translation sensor           handle SENSOR  
    end
    
    methods
        function FLOATJ = FLOATJOINT(Name,Type,Seg1,Seg2)
            FLOATJ = FLOATJ@JOINT(Name,Type,Seg1,Seg2); % Call JOINT constructor
        end
    end
    
end