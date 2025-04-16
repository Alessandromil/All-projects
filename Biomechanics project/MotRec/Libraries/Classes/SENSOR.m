classdef SENSOR < handle
    %SENSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name 
        Type
    end
    
    methods
        function S = SENSOR(Name,Type)
            S.Name = Name;
            S.Type = Type;
        end
    end
    
end

